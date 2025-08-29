#include "Primitive.hpp"
#include <chrono>
AABB GeometricPrimitive::Bounding_box() const {
    return shape->Bounding_box();
}
bool GeometricPrimitive::IntersectPred(const Ray& ray, float max) const {
    bool hit = shape->IntersectPred(ray,max);
    //if(material->Alpha()) we need UV coords from intersecpt pred!!
    return hit;
}
bool GeometricPrimitive::Intersect(const Ray& ray, SurfaceInteraction& interaction, float max) const {
    bool hit = shape->Intersect(ray,interaction,max);
    if(hit){
        interaction.AreaLight = areaLight;
        interaction.mat = material;
    }
    return hit;
}

std::vector<Light*> GeometricPrimitive::GetLights() const {
    return areaLight != nullptr ? std::vector<Light*>{areaLight} : std::vector<Light*>{};
}
/*
BVH::BVH(const std::vector<GeometricPrimitive>& prims){
    auto start = std::chrono::high_resolution_clock::now();
    nodes.reserve(prims.size()*2-1);
    primitives.reserve(prims.size());
    std::vector<PrimitiveInfo> primitiveInfo;
    primitiveInfo.reserve(prims.size());
    for(int i = 0;i<prims.size();i++){
        AABB bbox = prims[i].Bounding_box();
        primitiveInfo.emplace_back(i,bbox);
    }
    build_bvh(0,prims.size(),primitiveInfo);
    this->nodes.shrink_to_fit();
    for(const auto& info : primitiveInfo){
        primitives.emplace_back(prims[info.index]);
    }
    auto duration = std::chrono::high_resolution_clock::now() - start;
    std::cout<<"BVH build time in ms: "<<std::chrono::duration_cast<std::chrono::milliseconds>(duration).count()<<"\n";
}

bool BVH::Intersect(const Ray& ray, SurfaceInteraction& interaction, float max) const {
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


bool BVH::IntersectPred(const Ray& ray, float max) const {
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

inline bool BVH::intersectPrimitivesPred(const Ray& ray, float max, int right, int count) const{
    for(int triangle_index = right;triangle_index<right+count;triangle_index++){
        if(primitives[triangle_index].IntersectPred(ray,max))return true;
    }
    return false;
}

inline bool BVH::intersectPrimitives(const Ray& ray, float& max, int right, int count, SurfaceInteraction& interaction) const {
    bool hit = false;
    for(int triangle_index = right;triangle_index<right+count;triangle_index++){
        if(primitives[triangle_index].Intersect(ray,interaction,max)){
            hit = true;
            max = interaction.t;
        }
    }
    
    return hit;
}
*/