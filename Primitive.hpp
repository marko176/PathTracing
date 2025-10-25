#pragma once
#include "Shape.hpp"
#include "Light.hpp"
#include "Medium.hpp"
#include <chrono>

class Primitive {
public:
    virtual ~Primitive() = default;
    virtual AABB BoundingBox() const = 0;
    virtual bool IntersectPred(const Ray& ray, float max = std::numeric_limits<float>::infinity()) const = 0;
    virtual bool Intersect(const Ray& ray, SurfaceInteraction& interaction, float max = std::numeric_limits<float>::infinity()) const = 0;
    virtual std::vector<std::shared_ptr<Light>> GetLights() const = 0;
};

//chech with all shared ptr
class GeometricPrimitive : public Primitive{
public:
    virtual ~GeometricPrimitive() = default;
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
    virtual ~TransformedPrimitive() = default;
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



class AnimatedPrimitive : public Primitive {
public:
    virtual ~AnimatedPrimitive() = default;
    AnimatedPrimitive(const std::shared_ptr<Primitive>& primitive,const glm::vec3& direction,const glm::vec2& timeBounds) : primitive(primitive), dir(direction), timeBounds(timeBounds) {

    }
    AABB BoundingBox() const final ;
    bool IntersectPred(const Ray& ray, float max = std::numeric_limits<float>::infinity()) const final ;
    bool Intersect(const Ray& ray, SurfaceInteraction& interaction, float max = std::numeric_limits<float>::infinity()) const final ;
    std::vector<std::shared_ptr<Light>> GetLights() const final ;
private:
    std::shared_ptr<Primitive> primitive;
    glm::vec3 dir;
    glm::vec2 timeBounds;
};



struct Bin {
    AABB aabb;
    uint32_t triCount = 0;
};

struct alignas(32) BVH_NODE{
    AABB bbox;
    uint32_t right = 0;// right node / first triangle / meshID  / modelID
    uint32_t count = 0;// is_leaf / tirangle count / mesh count / modelCount
    //right holds first triangle index
};

struct alignas(32) simdBVH_NODE{
    alignas(16) float min[3];
    float right = 0;// right node / first triangle / meshID  / modelID
    alignas(16) float max[3];
    float count = 0;// is_leaf / tirangle count / mesh count / modelCount
    //right holds first triangle index
};



//maybe have BVHBase 

template <typename T>
class BVH : public Primitive{
public:
    virtual ~BVH() = default;
    BVH() = default;

    BVH(const std::vector<T>& prims) {
        auto start = std::chrono::high_resolution_clock::now();
        std::cout<<"Building BVH"<<std::endl;
        
        nodes.reserve(prims.size()*2-1);
        primitives.reserve(prims.size());

        std::vector<PrimitiveInfo> primitiveInfo;
        primitiveInfo.reserve(prims.size());

        for(std::size_t i = 0;i<prims.size();i++){
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
        std::cout<<"BVH built in: "<<std::chrono::duration_cast<std::chrono::milliseconds>(duration)<<std::endl;

        if(nodes.empty()){
            BVHbbox = AABB{};
        }else{
            BVHbbox = nodes[0].bbox;
        }
    }


    bool IntersectPred(const Ray& ray, float max = std::numeric_limits<float>::infinity()) const final {
        if(!BVHbbox.Hit(ray,max))return false;
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
        if(!BVHbbox.Hit(ray,max))return false;
        uint32_t stack[32];
        int i = 0;
        stack[i++]=0;
        
        bool hit_anything = false;
        
        while(i){
            int index = stack[--i];
            const BVH_NODE& node = nodes[index];
            
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
        return BVHbbox;
    }

    
private:

    struct PrimitiveInfo{
        PrimitiveInfo(uint32_t index, const AABB& bbox) : index{index}, bbox{bbox}, centroid(0.5f * (bbox.max + bbox.min)) {

        }

        uint32_t index;
        AABB bbox;
        glm::vec3 centroid;
    };

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
                bins.assign(BINS,Bin{});
                rightArea.clear();
                
                
                float max = -std::numeric_limits<float>::max();
                float min = std::numeric_limits<float>::max();
                for(int i = first_triangle;i<last_triangle;i++){
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
            
            int mid = std::partition(primitiveInfo.begin() + first_triangle,primitiveInfo.begin() + last_triangle,[&](const PrimitiveInfo& info){
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

    Mesh* mesh;
    AABB BVHbbox;
    std::vector<BVH_NODE> nodes;
    std::vector<T> primitives;
};

using BLAS = BVH<GeometricPrimitive>;
using TLAS = BVH<std::shared_ptr<Primitive>>;

template <typename T>
class simdBVH : public Primitive{
public:
    virtual ~simdBVH() = default;
    simdBVH() = default;

    simdBVH(const std::vector<T>& prims) {
        auto start = std::chrono::high_resolution_clock::now();
        std::cout<<"Building BVH"<<std::endl;
        nodes.reserve(prims.size()*2-1);
        primitives.reserve(prims.size());
        std::vector<PrimitiveInfo> primitiveInfo;
        primitiveInfo.reserve(prims.size());
        for(std::size_t i = 0;i<prims.size();i++){
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
        std::cout<<"BVH built in: "<<std::chrono::duration_cast<std::chrono::milliseconds>(duration)<<std::endl;
        if(nodes.empty()){
            BVHbbox = AABB{};
        }else{
            BVHbbox = AABB{{nodes[0].min[0],nodes[0].min[1],nodes[0].min[2]}};
            BVHbbox.Expand({nodes[0].max[0],nodes[0].max[1],nodes[0].max[2]});
        }
    }


    static inline float IntersectAABB(const Ray& ray,const float* min4,const float* max4, float max_t) {
        const __m128 bmin4 = _mm_load_ps(min4);
        const __m128 bmax4 = _mm_load_ps(max4);

        const __m128 t1 = _mm_mul_ps(_mm_sub_ps(bmin4,ray.O4),ray.rD4 );
        const __m128 t2 = _mm_mul_ps(_mm_sub_ps(bmax4,ray.O4),ray.rD4 );
        const __m128 vmax4 = _mm_max_ps(t1,t2);
        const __m128 vmin4 = _mm_min_ps(t1,t2);

        alignas(16) float temp_min[4],temp_max[4];
        _mm_store_ps(temp_min,vmin4);
        _mm_store_ps(temp_max,vmax4);

        const float tEntry = std::max(std::max(temp_min[0], temp_min[1]), temp_min[2]);
        const float tExit  = std::min(std::min(temp_max[0], temp_max[1]), temp_max[2]);

        return (tEntry <= tExit && tExit >= shadowEpsilon && tEntry <= max_t ) ? tEntry : std::numeric_limits<float>::infinity();//maybe texit greater than 0 ?
    }

    static inline bool HitAABB(const Ray& ray,const float* min4,const float* max4, float max_t) {
        return IntersectAABB(ray,min4,max4,max_t) != std::numeric_limits<float>::infinity();//maybe texit greater than 0 ?
    }


    bool IntersectPred(const Ray& ray, float max = std::numeric_limits<float>::infinity()) const final {
        if(!HitAABB(ray,nodes[0].min,nodes[0].max,max))return false;
        uint32_t stack[32];
        int i = 0;
        stack[i++]=0;
        
        while(i){
            int index = stack[--i];
            const simdBVH_NODE& node = nodes[index];
            
            // switch to intercesion test
            
            if(node.count == 0){
                int child1 = index+1;
                int child2 = node.right;
                if(HitAABB(ray,nodes[child1].min,nodes[child1].max,max))stack[i++]=child1;
                if(HitAABB(ray,nodes[child2].min,nodes[child2].max,max))stack[i++]=child2;
            }else if(intersectPrimitivesPred(ray,max,node.right,node.count)){
                return true;
            }
            
        }
        
        return false;
    }


    bool Intersect(const Ray& ray, SurfaceInteraction& interaction, float max = std::numeric_limits<float>::infinity()) const final {
        //will fail on empty nodes!
        if(!HitAABB(ray,nodes[0].min,nodes[0].max,max))return false;
        uint32_t stack[32];
        int i = 0;
        stack[i++]=0;
        
        bool hit_anything = false;
        
        while(i){
            int index = stack[--i];
            const simdBVH_NODE& node = nodes[index];
            
            if(node.count == 0){
                int child1 = index+1;
                int child2 = node.right;
                float dist_1 = IntersectAABB(ray,nodes[child1].min,nodes[child1].max,max);
                float dist_2 = IntersectAABB(ray,nodes[child2].min,nodes[child2].max,max);
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
        return BVHbbox;
    }

    
private:

    struct PrimitiveInfo{
        PrimitiveInfo(uint32_t index, const AABB& bbox) : index{index}, bbox{bbox}, centroid(0.5f * (bbox.max + bbox.min)) {

        }

        uint32_t index;
        AABB bbox;
        glm::vec3 centroid;
    };

    int build_bvh(int first_triangle, int last_triangle,std::vector<PrimitiveInfo>& primitiveInfo){
        int index = nodes.size();
        simdBVH_NODE& node = nodes.emplace_back();
        size_t object_span = last_triangle - first_triangle;

        AABB bbox;
        for(int i = first_triangle;i<last_triangle;i++){
            bbox.Expand(primitiveInfo[i].bbox);
        }
        node.min[0]=bbox.min[0];
        node.min[1]=bbox.min[1];
        node.min[2]=bbox.min[2];
        node.max[0]=bbox.max[0];
        node.max[1]=bbox.max[1];
        node.max[2]=bbox.max[2];
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
                bins.assign(BINS,Bin{});
                rightArea.clear();
                
                
                float max = -std::numeric_limits<float>::max();
                float min = std::numeric_limits<float>::max();
                for(int i = first_triangle;i<last_triangle;i++){
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
            float parent_cost = bbox.Area() * object_span;
            if(bestCost >= parent_cost){
                return index;
            }
            
            int mid = std::partition(primitiveInfo.begin() + first_triangle,primitiveInfo.begin() + last_triangle,[&](const PrimitiveInfo& info){
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

    Mesh* mesh;
    AABB BVHbbox;
    std::vector<simdBVH_NODE> nodes;
    std::vector<T> primitives;
};

using simdBLAS = simdBVH<GeometricPrimitive>;
using simdTLAS = simdBVH<std::shared_ptr<Primitive>>;