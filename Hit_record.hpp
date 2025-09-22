#pragma once
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtx/intersect.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include "Ray.hpp"
#include <memory>

class Material;
class Medium;
struct hittable;
struct LightObject;

class Light;

struct GeometricInteraction {
    glm::vec3 p = {0,0,0};
    glm::vec3 n = {0,0,0};//geometric normal
    std::shared_ptr<Medium> medium = nullptr;
    
    bool isSurfaceInteraction() const {
        return n != glm::vec3(0,0,0);
    }

    bool isMediumInteraction() const {
        return !isSurfaceInteraction();
    }

    std::shared_ptr<Medium> getMedium(const glm::vec3& dir) const{
        if(glm::dot(dir,n)<0)return medium;//going into surface
        return nullptr;
    }
    GeometricInteraction() = default;
    GeometricInteraction(const glm::vec3& point, const glm::vec3& geometric_normal,const std::shared_ptr<Medium>& medium = nullptr) : p{point} , n{geometric_normal} , medium(medium) {

    }
};

struct SurfaceInteraction : public GeometricInteraction{
    SurfaceInteraction(const glm::vec3& point, const glm::vec3& geometric_normal, const glm::vec2& uv,const std::shared_ptr<Medium>& medium = nullptr) : GeometricInteraction{point,geometric_normal,medium}, uv{uv} {

    }

    SurfaceInteraction() = default;
    SurfaceInteraction(const GeometricInteraction& interaction) : GeometricInteraction{interaction} {

    }
    glm::vec2 uv = {0,0};
    float t;//we dont need t if we have point? t = glm::length(ray)
    glm::vec3 ns;
    std::shared_ptr<Material> mat; // switch to material* no need for shared
    glm::vec3 tangent;
    glm::vec3 bitangent;
    std::shared_ptr<Light> AreaLight;
};

struct MediumInteraction : public GeometricInteraction {
    MediumInteraction(const glm::vec3& point, const glm::vec3& geometric_normal, const std::shared_ptr<Medium>& medium = nullptr) : GeometricInteraction{point,geometric_normal,medium} {

    }

    MediumInteraction() = default;

    bool isValid() const {
        //phasefunc != null
        return valid;
    }
    bool valid = false;
    //phase function
};


struct AABB{
    glm::vec3 min;
    glm::vec3 max;

    AABB() : min{std::numeric_limits<float>::max()} , max{-min} {

    }


    AABB(const glm::vec3& point) : min(point) , max(point) {

    }

    inline void expand(const AABB& other) {
        min = glm::min(min,other.min);
        max = glm::max(max,other.max);
    }

    inline void expand(const glm::vec3& point) {
        min = glm::min(min,point);
        max = glm::max(max,point);
    }

    inline bool hit(const Ray& ray,float min_t,float max_t) const {// add interval
        glm::vec3 invD = ray.inv_dir;//store in ray!!!! test speed
        // Calculate t for each slab
        glm::vec3 t0s = (min - ray.origin) * invD;
        glm::vec3 t1s = (max - ray.origin) * invD;
        // Component‐wise min/max
        glm::vec3 tSmaller = glm::min(t0s, t1s);
        glm::vec3 tLarger  = glm::max(t0s, t1s);
        // Over all axes
        float tEntry = std::max(std::max(tSmaller.x, tSmaller.y), tSmaller.z);
        float tExit  = std::min(std::min(tLarger.x,  tLarger.y),  tLarger.z);
        return (tEntry <= tExit) && (tExit >= min_t) && (tEntry <= max_t);
    }

    inline bool hit(const Ray& ray,float max_t) const {// add interval
        glm::vec3 invD = ray.inv_dir;//store in ray!!!! test speed
        // Calculate t for each slab
        glm::vec3 t0s = (min - ray.origin) * invD;
        glm::vec3 t1s = (max - ray.origin) * invD;
        // Component‐wise min/max
        glm::vec3 tSmaller = glm::min(t0s, t1s);
        glm::vec3 tLarger  = glm::max(t0s, t1s);
        // Over all axes
        float tEntry = std::max(std::max(tSmaller.x, tSmaller.y), tSmaller.z);
        float tExit  = std::min(std::min(tLarger.x,  tLarger.y),  tLarger.z);
        return (tEntry <= tExit) && (tExit > 0.001f) && (tEntry <= max_t);
    }

    inline float hit_other(const Ray& ray,float min_t,float max_t) const {// add interval
        glm::vec3 invD = ray.inv_dir;//store in ray!!!! test speed
        // Calculate t for each slab
        glm::vec3 t0s = (min - ray.origin) * invD;
        glm::vec3 t1s = (max - ray.origin) * invD;
        // Component‐wise min/max
        glm::vec3 tSmaller = glm::min(t0s, t1s);
        glm::vec3 tLarger  = glm::max(t0s, t1s);
        // Over all axes
        float tEntry = std::max(std::max(tSmaller.x, tSmaller.y), tSmaller.z);
        float tExit  = std::min(std::min(tLarger.x,  tLarger.y),  tLarger.z);
        return (tEntry <= tExit && tEntry <= max_t && tExit >= min_t) ? tEntry : std::numeric_limits<float>::infinity();
    }

    inline float hit_other(const Ray& ray,float max_t) const {// add interval
        glm::vec3 invD = ray.inv_dir;//store in ray!!!! test speed
        // Calculate t for each slab
        glm::vec3 t0s = (min - ray.origin) * invD;
        glm::vec3 t1s = (max - ray.origin) * invD;
        // Component‐wise min/max
        glm::vec3 tSmaller = glm::min(t0s, t1s);
        glm::vec3 tLarger  = glm::max(t0s, t1s);
        // Over all axes
        float tEntry = std::max(std::max(tSmaller.x, tSmaller.y), tSmaller.z);
        float tExit  = std::min(std::min(tLarger.x,  tLarger.y),  tLarger.z);
        return (tEntry <= tExit && tEntry <= max_t && tExit > 0.001f) ? tEntry : std::numeric_limits<float>::infinity();
    }

    

    inline float area() const { 
        glm::vec3 e = max - min; // box extent
        return e.x * e.y + e.y * e.z + e.z * e.x; 
    }

    glm::vec3 operator[](int i) const { return (i == 0) ? min : max; }
    glm::vec3 &operator[](int i) { return (i == 0) ? min : max; }

    glm::vec3 Corner(int n){
        return {(*this)[(n & 1)].x,
                (*this)[(n & 2) ? 1 : 0].y,
                (*this)[(n & 4) ? 1 : 0].z};
    }
};

struct Bounds2i {
    glm::ivec2 min = {0,0};
    glm::ivec2 max = {0,0};
};








