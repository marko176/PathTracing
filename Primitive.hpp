#pragma once
#include "Shape.hpp"
#include "Light.hpp"
#include <chrono>
class Primitive {
public:
    virtual AABB Bounding_box() const = 0;
    virtual bool IntersectPred(const Ray& ray, float max = 1e30f) const = 0;
    virtual bool Intersect(const Ray& ray, SurfaceInteraction& interaction, float max = 1e30f) const = 0;
    virtual std::vector<Light*> GetLights() const = 0;
};

class GeometricPrimitive : public Primitive{
public:

    GeometricPrimitive(Shape* primitive_shape, const std::shared_ptr<Material>& material,AreaLight* areaLight) : shape{primitive_shape} , material{material} , areaLight{areaLight} {

    }

    virtual AABB Bounding_box() const final ;
    virtual bool IntersectPred(const Ray& ray, float max = 1e30f) const final ;
    virtual bool Intersect(const Ray& ray, SurfaceInteraction& interaction, float max = 1e30f) const final ;
    virtual std::vector<Light*> GetLights() const final ;
private:

    Shape* shape;
    std::shared_ptr<Material> material;
    AreaLight* areaLight;
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

    Mesh* mesh;
    std::vector<BVH_NODE> nodes;
    std::vector<T> primitives;

    BVH() = default;

    BVH(const std::vector<T>& prims) {
        auto start = std::chrono::high_resolution_clock::now();
        nodes.reserve(prims.size()*2-1);
        primitives.reserve(prims.size());
        std::vector<PrimitiveInfo> primitiveInfo;
        primitiveInfo.reserve(prims.size());
        for(int i = 0;i<prims.size();i++){
            if constexpr(std::is_pointer_v<T>){
                AABB bbox = prims[i]->Bounding_box();
                primitiveInfo.emplace_back(i,bbox);
            }else{
                AABB bbox = prims[i].Bounding_box();
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
            node.bbox.expand(primitiveInfo[i].bbox);
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
                    bins[binId].aabb.expand(primitiveInfo[i].bbox);

                }

                AABB leftBox,rightBox;
                for(int i = 0;i<BINS-1;i++){
                    rightBox.expand( bins[BINS - 1 - i].aabb );
                    rightArea[BINS - 2 - i] = rightBox.area();
                }

                scale = (max-min)/BINS;
                int leftSum = 0;
                for(int i = 0;i<BINS-1;i++){
                    leftSum += bins[i].triCount;
                    leftBox.expand( bins[i].aabb );
                    float cost = leftSum * leftBox.area() + (object_span - leftSum) * rightArea[i];
                    if(cost < bestCost){
                        best_axis = axis;
                        bestPos = min + (i+1) * scale;
                        bestCost = cost;
                    }
                }

            }
            float parent_cost = node.bbox.area() * object_span;
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

    bool IntersectPred(const Ray& ray, float max = 1e30f) const final {
        if(!nodes[0].bbox.hit(ray,max))return false;
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
                if(nodes[child1].bbox.hit(ray,max))stack[i++]=child1;
                if(nodes[child2].bbox.hit(ray,max))stack[i++]=child2;
            }else if(intersectPrimitivesPred(ray,max,node.right,node.count)){
                return true;
            }
            
        }
        
        return false;
    }


    bool Intersect(const Ray& ray, SurfaceInteraction& interaction, float max = 1e30f) const final {
        if(!nodes[0].bbox.hit(ray,max))return false;
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
                float dist_1 = nodes[child1].bbox.hit_other(ray,max);
                float dist_2 = nodes[child2].bbox.hit_other(ray,max);
                if(dist_1 > dist_2){
                    std::swap(dist_1,dist_2);
                    std::swap(child1,child2);
                }
                if(dist_2 != 1e30f){
                    stack[i++]=child2;
                    stack[i++]=child1;
                }else if(dist_1 != 1e30f){
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
            if constexpr(std::is_pointer_v<T>){
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
            if constexpr(std::is_pointer_v<T>){
                if(primitives[triangle_index]->IntersectPred(ray,max))return true;
            }else {
                if(primitives[triangle_index].IntersectPred(ray,max))return true;
            }
        }
        return false;
    }

    
    std::vector<Light*> GetLights() const final {
        std::vector<Light*> lights;
        for(const T& prim : primitives){
            std::vector<Light*> primLights;
            if constexpr(std::is_pointer_v<T>){
                primLights = prim->GetLights();
            }else {
                primLights = prim.GetLights();
            }
            lights.insert(lights.end(),primLights.begin(),primLights.end());
        }
        return lights;
    }

    AABB Bounding_box() const final{
        return nodes.empty() ? AABB{} : nodes[0].bbox;
    }

    //store nodes in flat array
    //nodes hole left index, right index, is_leaf and ptr to hittable -> same as before but flat
   
    
};

using BLAS = BVH<GeometricPrimitive>;
using TLAS = BVH<Primitive*>;

struct TLAS_BVH_Prim : public Primitive{


    struct Bin {
        AABB aabb;
        uint32_t triCount = 0;
    };

    std::vector<Primitive*> hittables;//blas
    std::vector<BVH_NODE> nodes;
    TLAS_BVH_Prim() = default;//should switch to method called setup()
    TLAS_BVH_Prim(const std::vector<Primitive*>& hittables) : hittables(hittables) {
        nodes.reserve(hittables.size()*2-1);
        build_bvh(0,hittables.size());
        this->nodes.shrink_to_fit();
        std::cout<<"TLAS BUILT\n";
    }

    bool add(Primitive* prim) {
        hittables.push_back(prim);
        return true;
    }

    bool build() {
        nodes.reserve(hittables.size()*2-1);
        build_bvh(0,hittables.size());
        this->nodes.shrink_to_fit();
        std::cout<<"TLAS BUILT\n";
        return true;
    }



    int build_bvh(int first_triangle, int last_triangle){
        int index = nodes.size();
        BVH_NODE& node = nodes.emplace_back();
        for(int i = first_triangle;i<last_triangle;i++){
            node.bbox.expand(hittables[i]->Bounding_box());
        }

      
        


        size_t object_span = last_triangle - first_triangle;
        node.count = 0;
        if (object_span <= 2) {
            node.right = first_triangle;//triangle
            node.count = object_span;
        } else {
            auto evaluateSAH = [&](int first, int last, int axis,float pos){
                AABB left,right;
                int left_count = 0;
                int right_count = 0;
                for(int i = first;i<last;i++){
                    float triangle_centroid = (hittables[i]->Bounding_box().max[axis]+hittables[i]->Bounding_box().min[axis])/2;
                    if(triangle_centroid < pos){
                        left_count++;
                        left.expand(hittables[i]->Bounding_box());
                    }else{
                        right_count++;
                        right.expand(hittables[i]->Bounding_box());
                    }
                }
                float cost = left_count * left.area() + right_count * right.area();
                return cost > 0 ? cost : 1e30f;
            };
            int best_axis = 0;
            float bestPos = 0;
            float bestCost = 1e30f;
            for(int axis = 0;axis<3;axis++){
          
                    int count = 100;
                    float scale = (node.bbox.max[axis] - node.bbox.min[axis])/count;
                    for(int i = 0;i<count;i++){
                        
                        float cost = evaluateSAH(first_triangle,last_triangle,axis,node.bbox.min[axis] + i * scale);
                        if(cost < bestCost){
                            bestPos = node.bbox.min[axis] + i * scale;
                            best_axis = axis;
                            bestCost = cost;
                        }
                    }
                
        
            }
       

            
            int mid = std::partition(hittables.begin() + first_triangle,hittables.begin() + last_triangle,[&](Primitive* h){
                 float triangle_centroid = (h->Bounding_box().max[best_axis]+h->Bounding_box().min[best_axis])/2;
                return triangle_centroid < bestPos;
            }) - hittables.begin();
        
            //node.left = pos+1;
            build_bvh(first_triangle, mid);
            node.right = build_bvh(mid, last_triangle);
        }
        return index;
    }
    AABB Bounding_box() const override{
        return nodes[0].bbox;
    }

    std::vector<Light*> GetLights() const override {
        std::vector<Light*> lights;
        for(auto* prim : hittables){
            if(prim == nullptr)continue;
            std::vector<Light*> primLights = prim->GetLights();
            lights.insert(lights.end(),primLights.begin(),primLights.end());
        }
        return lights;
    }

    virtual bool Intersect(const Ray& ray, SurfaceInteraction& interaction, float max = 1e30f) const override {
        if(!nodes[0].bbox.hit(ray,max))return false;
        uint32_t stack[24];
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
                float dist_1 = nodes[child1].bbox.hit_other(ray,max);
                float dist_2 = nodes[child2].bbox.hit_other(ray,max);
                if(dist_1 > dist_2){
                    std::swap(dist_1,dist_2);
                    std::swap(child1,child2);
                }
                if(dist_2 != 1e30f){
                    stack[i++]=child2;
                    stack[i++]=child1;
                }else if(dist_1 != 1e30f){
                    stack[i++]=child1;
                }
         
            }else if(hit_hittable(ray,max,node.right,node.count,interaction)){
                hit_anything = true;
            }
            
        }
        
        return hit_anything;
    }

    virtual bool IntersectPred(const Ray& ray, float max = 1e30f) const override{
        if(!nodes[0].bbox.hit(ray,max))return false;
        uint32_t stack[24];
        int i = 0;
        stack[i++]=0;
        
        
        while(i){
            int index = stack[--i];
            const BVH_NODE& node = nodes[index];
            
            // switch to intercesion test
            
            if(node.count == 0){
                int child1 = index+1;
                int child2 = node.right;
                if(nodes[child1].bbox.hit(ray,max))stack[i++]=child1;
                if(nodes[child2].bbox.hit(ray,max))stack[i++]=child2;
         
            }else if(hit_hittablePred(ray,max,node.right,node.count)){
                return true;
            }
            
        }
        
        return false;
    }


    inline bool hit_hittable(const Ray& ray, float& max,int right, int count, SurfaceInteraction& interaction) const{
        bool hit = false;
        for(int i = right;i<right+count;i++){
            if(hittables[i]->Intersect(ray,interaction,max)){
                hit = true;
                max = interaction.t;
            }
        }
        return hit;
    }  

    inline bool hit_hittablePred(const Ray& ray,float max, int right, int count) const{
        for(int i = right;i<right+count;i++){
            if(hittables[i]->IntersectPred(ray,max))
                return true;
        }
        return false;
    }  
    
 
};