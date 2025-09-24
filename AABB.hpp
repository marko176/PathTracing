#pragma once
#include "Ray.hpp"
#include "Interaction.hpp"





struct AABB{
    glm::vec3 min;
    glm::vec3 max;

    AABB() : min{std::numeric_limits<float>::max()} , max{-min} {

    }


    AABB(const glm::vec3& point) : min(point) , max(point) {

    }

    inline void Expand(const AABB& other) {
        min = glm::min(min,other.min);
        max = glm::max(max,other.max);
    }

    inline void Expand(const glm::vec3& point) {
        min = glm::min(min,point);
        max = glm::max(max,point);
    }

    inline bool Hit(const Ray& ray,float min_t,float max_t) const {// add interval
        return HitDistance(ray,min_t,max_t) != std::numeric_limits<float>::infinity();
    }

    inline bool Hit(const Ray& ray,float max_t) const {// add interval
        return HitDistance(ray,max_t) != std::numeric_limits<float>::infinity();
    }

    inline float HitDistance(const Ray& ray,float min_t,float max_t) const {// add interval
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

    inline float HitDistance(const Ray& ray,float max_t) const {// add interval
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
        return (tEntry <= tExit && tEntry <= max_t && tExit > 0) ? tEntry : std::numeric_limits<float>::infinity();//maybe texit greater than 0 ?
    }

   
    //using AABB::nothing
    //AABB::universe
    //mediumintr takes phasefunc

    inline float Area() const { 
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








