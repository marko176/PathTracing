#pragma once
#include "Shape.hpp"
#include "Light.hpp"
#include "Medium.hpp"
#include <chrono>

class Primitive {
public:
    virtual AABB BoundingBox() const = 0;
    virtual bool IntersectPred(const Ray& ray, float max = std::numeric_limits<float>::infinity()) const = 0;
    virtual bool Intersect(const Ray& ray, SurfaceInteraction& interaction, float max = std::numeric_limits<float>::infinity()) const = 0;
    virtual std::vector<std::shared_ptr<Light>> GetLights() const = 0;
};

//chech with all shared ptr
class GeometricPrimitive : public Primitive{
public:
 
    GeometricPrimitive(const std::shared_ptr<Shape>& primitive_shape, const std::shared_ptr<Material>& material,const std::shared_ptr<AreaLight>& areaLight = nullptr,const std::shared_ptr<Medium>& medium = nullptr) : shape{primitive_shape} , material{material} , areaLight{areaLight} , medium(medium) {}
    
    AABB BoundingBox() const final ;
    bool IntersectPred(const Ray& ray, float max = std::numeric_limits<float>::infinity()) const final ;
    bool Intersect(const Ray& ray, SurfaceInteraction& interaction, float max = std::numeric_limits<float>::infinity()) const final ;
    std::vector<std::shared_ptr<Light>> GetLights() const final ;
private:
    std::shared_ptr<Shape> shape;
    std::shared_ptr<Material> material;
    std::shared_ptr<AreaLight> areaLight;
    std::shared_ptr<Medium> medium;
};


class TransformedPrimitive : public Primitive {
public:
    TransformedPrimitive(const std::shared_ptr<Primitive>& primitive,const glm::mat4& transform) : primitive(primitive), transform(transform) , invTransform(glm::inverse(transform)) {

    }
    AABB BoundingBox() const final ;
    bool IntersectPred(const Ray& ray, float max = std::numeric_limits<float>::infinity()) const final ;
    bool Intersect(const Ray& ray, SurfaceInteraction& interaction, float max = std::numeric_limits<float>::infinity()) const final ;
    std::vector<std::shared_ptr<Light>> GetLights() const final ;
private:
    std::shared_ptr<Primitive> primitive;
    glm::mat4 transform;
    glm::mat4 invTransform;
};

struct BVH_NODE{
    AABB bbox;
    uint32_t right = 0;// right node / first triangle / meshID  / modelID
    uint32_t count = 0;// is_leaf / tirangle count / mesh count / modelCount
    //right holds first triangle index
};


template <typename T>
struct BVH : public Primitive{

    struct Bin {
        AABB aabb;
        uint32_t triCount = 0;
    };

    struct PrimitiveInfo{
        PrimitiveInfo(uint32_t index, const AABB& bbox) : index{index}, bbox{bbox}, centroid(0.5f * (bbox.max + bbox.min)) {

        }

        uint32_t index;
        AABB bbox;
        glm::vec3 centroid;
    };

    
    
    BVH() = default;

    BVH(const std::vector<T>& prims) {
        auto start = std::chrono::high_resolution_clock::now();
        nodes.reserve(prims.size()*2-1);
        primitives.reserve(prims.size());
        std::vector<PrimitiveInfo> primitiveInfo;
        primitiveInfo.reserve(prims.size());
        for(int i = 0;i<prims.size();i++){
            if constexpr(requires { prims[i]->BoundingBox(); }){
                AABB bbox = prims[i]->BoundingBox();
                primitiveInfo.emplace_back(i,bbox);
            }else{
                AABB bbox = prims[i].BoundingBox();
                primitiveInfo.emplace_back(i,bbox);
            }
        }
        build_bvh(0,prims.size(),primitiveInfo);
        this->nodes.shrink_to_fit();
        for(const auto& info : primitiveInfo){
            primitives.emplace_back(prims[info.index]);
        }
        auto duration = std::chrono::high_resolution_clock::now() - start;
        std::cout<<"BVH build time in ms: "<<std::chrono::duration_cast<std::chrono::milliseconds>(duration).count()<<"\n";
    }

    
    int build_bvh(int first_triangle, int last_triangle,std::vector<PrimitiveInfo>& primitiveInfo){
        int index = nodes.size();
        BVH_NODE& node = nodes.emplace_back();
        size_t object_span = last_triangle - first_triangle;
        for(int i = first_triangle;i<last_triangle;i++){
            node.bbox.Expand(primitiveInfo[i].bbox);
        }



       
        node.right = first_triangle;//triangle
        node.count = object_span;
        if (object_span > 2){

            int best_axis = 0;
            float bestPos = 0;
            float bestCost = std::numeric_limits<float>::max();
            int BINS =  object_span >= 8192 ? 64 :
                        object_span >= 1024 ? 32 :
                        object_span >= 64 ? 16 : 8;
            std::vector<Bin> bins;
            std::vector<float> rightArea;
            bins.reserve(BINS);
            rightArea.reserve(BINS-1);
            for(int axis = 0;axis<3;axis++){
                //int count = 200;//set count to sqrt(span)
                bins.assign(BINS,Bin{});
                rightArea.clear();
                
                float max = -std::numeric_limits<float>::max();
                float min = std::numeric_limits<float>::max();
                for(int i = first_triangle;i<last_triangle;i++){
                    //float triangle_centroid = (vertices[indices[i*3+0]][axis]+vertices[indices[i*3+1]][axis]+vertices[indices[i*3+2]][axis])/3.0f;
                    float triangle_centroid = primitiveInfo[i].centroid[axis];
                    max = std::max(max,triangle_centroid);
                    min = std::min(min,triangle_centroid);
                }
                if(max == min)continue;
                float scale = BINS/(max - min);

                for(int i = first_triangle;i<last_triangle;i++){
                    //float triangle_centroid = (vertices[indices[i*3+0]][axis]+vertices[indices[i*3+1]][axis]+vertices[indices[i*3+2]][axis])/3.0f;
                    float triangle_centroid = primitiveInfo[i].centroid[axis];
                    int binId = std::min<int>(BINS-1, (triangle_centroid - min)*scale);
                    ++bins[binId].triCount;
                    bins[binId].aabb.Expand(primitiveInfo[i].bbox);

                }

                AABB leftBox,rightBox;
                for(int i = 0;i<BINS-1;i++){
                    rightBox.Expand( bins[BINS - 1 - i].aabb );
                    rightArea[BINS - 2 - i] = rightBox.Area();
                }

                scale = (max-min)/BINS;
                int leftSum = 0;
                for(int i = 0;i<BINS-1;i++){
                    leftSum += bins[i].triCount;
                    leftBox.Expand( bins[i].aabb );
                    float cost = leftSum * leftBox.Area() + (object_span - leftSum) * rightArea[i];
                    if(cost < bestCost){
                        best_axis = axis;
                        bestPos = min + (i+1) * scale;
                        bestCost = cost;
                    }
                }

            }
            float parent_cost = node.bbox.Area() * object_span;
            if(bestCost >= parent_cost){
                return index;
            }

            int mid = std::partition(primitiveInfo.begin() + first_triangle,primitiveInfo.begin() + last_triangle,[&](PrimitiveInfo& info){
                float triangle_centroid = info.centroid[best_axis];
                return triangle_centroid <= bestPos;
            }) - primitiveInfo.begin();


            if(mid == first_triangle || mid == last_triangle){
                return index;
            }


        
            //node.left = pos+1;
            build_bvh(first_triangle, mid, primitiveInfo);
            node.right = build_bvh(mid, last_triangle, primitiveInfo);
            node.count = 0;
        }
        return index;
    }

    bool IntersectPred(const Ray& ray, float max = std::numeric_limits<float>::infinity()) const final {
        if(!nodes[0].bbox.Hit(ray,max))return false;
        uint32_t stack[32];
        int i = 0;
        stack[i++]=0;
        
        while(i){
            int index = stack[--i];
            const BVH_NODE& node = nodes[index];
            
            // switch to intercesion test
            
            if(node.count == 0){
                int child1 = index+1;
                int child2 = node.right;
                if(nodes[child1].bbox.Hit(ray,max))stack[i++]=child1;
                if(nodes[child2].bbox.Hit(ray,max))stack[i++]=child2;
            }else if(intersectPrimitivesPred(ray,max,node.right,node.count)){
                return true;
            }
            
        }
        
        return false;
    }


    bool Intersect(const Ray& ray, SurfaceInteraction& interaction, float max = std::numeric_limits<float>::infinity()) const final {
        if(!nodes[0].bbox.Hit(ray,max))return false;
        uint32_t stack[32];
        int i = 0;
        stack[i++]=0;
        
        bool hit_anything = false;
        
        while(i){
            int index = stack[--i];
            const BVH_NODE& node = nodes[index];
            
            // switch to intercesion test
            
            if(node.count == 0){
                int child1 = index+1;
                int child2 = node.right;
                float dist_1 = nodes[child1].bbox.HitDistance(ray,max);
                float dist_2 = nodes[child2].bbox.HitDistance(ray,max);
                if(dist_1 > dist_2){
                    std::swap(dist_1,dist_2);
                    std::swap(child1,child2);
                }
                if(dist_2 != std::numeric_limits<float>::infinity()){
                    stack[i++]=child2;
                    stack[i++]=child1;
                }else if(dist_1 != std::numeric_limits<float>::infinity()){
                    stack[i++]=child1;
                }
        
            }else if(intersectPrimitives(ray,max,node.right,node.count,interaction)){
                hit_anything = true;
            }
            
        }
        
        return hit_anything;
    }



    inline bool intersectPrimitives(const Ray& ray, float& max, int right, int count, SurfaceInteraction& interaction) const {
        bool hit = false;
        for(int triangle_index = right;triangle_index<right+count;triangle_index++){
            if constexpr(requires { primitives[triangle_index]->Intersect(ray,interaction,max); }){
                if(primitives[triangle_index]->Intersect(ray,interaction,max)){
                    hit = true;
                    max = interaction.t;
                }
            }else {
                if(primitives[triangle_index].Intersect(ray,interaction,max)){
                    hit = true;
                    max = interaction.t;
                }
            }
        }
        
        return hit;
    }

    inline bool intersectPrimitivesPred(const Ray& ray, float max, int right, int count) const {
        for(int triangle_index = right;triangle_index<right+count;triangle_index++){
            if constexpr (requires { primitives[triangle_index]->IntersectPred(ray,max); }){
                if(primitives[triangle_index]->IntersectPred(ray,max))return true;
            }else {
                if(primitives[triangle_index].IntersectPred(ray,max))return true;
            }
        }
        return false;
    }

    
    std::vector<std::shared_ptr<Light>> GetLights() const final {
        std::vector<std::shared_ptr<Light>> lights;
        for(const T& prim : primitives){
            std::vector<std::shared_ptr<Light>> primLights;
            if constexpr (requires { prim->GetLights(); }){
                primLights = prim->GetLights();
            }else {
                primLights = prim.GetLights();
            }
            lights.insert(lights.end(),primLights.begin(),primLights.end());
        }
        return lights;
    }

    AABB BoundingBox() const final{
        return nodes.empty() ? AABB{} : nodes[0].bbox;
    }

private:
    Mesh* mesh;
    std::vector<BVH_NODE> nodes;
    std::vector<T> primitives;
    
};

using BLAS = BVH<GeometricPrimitive>;
using TLAS = BVH<std::shared_ptr<Primitive>>;
