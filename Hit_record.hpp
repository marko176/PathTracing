#pragma once
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtx/intersect.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include "Ray.hpp"
#include <memory>

class Material;

struct hittable;
struct LightObject;

class Light;

struct GeometricInteraction {
    glm::vec3 p = {0,0,0};
    glm::vec3 n = {0,0,0};//geometric normal
    glm::vec2 uv = {0,0};

    GeometricInteraction() = default;
    GeometricInteraction(const glm::vec3& point, const glm::vec3& geometric_normal, const glm::vec2& uv) : p{point} , n{geometric_normal} , uv{uv} {

    }
};

struct SurfaceInteraction : public GeometricInteraction{
    SurfaceInteraction(const glm::vec3& point, const glm::vec3& geometric_normal, const glm::vec2& uv) : GeometricInteraction{point,geometric_normal,uv} {

    }

    SurfaceInteraction() = default;
    SurfaceInteraction(const GeometricInteraction& interaction) : GeometricInteraction{interaction} {

    }

    float t;//we dont need t if we have point? t = glm::length(ray)
    glm::vec3 ns;
    std::shared_ptr<Material> mat; // switch to material* no need for shared
    glm::vec3 tangent;
    glm::vec3 bitangent;
    const Light* AreaLight;
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
        return (tEntry <= tExit && tEntry <= max_t && tExit >= min_t) ? tEntry : 1e30f;
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
        return (tEntry <= tExit && tEntry <= max_t && tExit > 0.001f) ? tEntry : 1e30f;
    }

    

    float min_(int n) const {
        if (n == 1) return min.y;
        if (n == 2) return min.z;
        return min.x;
    }

    inline float area() const { 
        glm::vec3 e = max - min; // box extent
        return e.x * e.y + e.y * e.z + e.z * e.x; 
    }
};

struct hittable {
    virtual ~hittable() = default;
    virtual bool hit(const Ray& ray, float min, float max, SurfaceInteraction& interaction) const = 0;
    virtual AABB bounding_box() const = 0;
    virtual bool is_bvh_node() const {
        return false;
    }
    virtual bool hit_bvh_info(const Ray& ray, float min, float max, SurfaceInteraction& interaction,int& bvh_tests,int& triangle_tests) const {
        return false;
    }

    virtual bool intersectPred(const Ray& ray, float min, float max, SurfaceInteraction& interaction) const {
        return hit(ray,min,max,interaction);
    }

    virtual float area() const {
        return 0;
    }

    virtual glm::vec3 sampleSurface() const {
        //should also return pdf if prob is not uniform
        return {0,0,0};
    }

    virtual glm::vec3 sampleSurface(const glm::vec2& UV) const {
        //should also return pdf if prob is not uniform
        return {0,0,0};
    }
};






