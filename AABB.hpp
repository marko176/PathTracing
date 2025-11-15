#pragma once
#include "Ray.hpp"
#include "Interaction.hpp"


constexpr inline float shadowEpsilon = 0.00001f;


struct AABB{
    glm::vec3 min;
    glm::vec3 max;



    AABB() : min { std::numeric_limits<float>::infinity() }, max { -min }{

    }


    AABB(const glm::vec3& point) : min(point), max(point){

    }

    inline void Expand(const AABB& other){
        min = glm::min(min, other.min);
        max = glm::max(max, other.max);
    }

    inline void Expand(const glm::vec3& point){
        min = glm::min(min, point);
        max = glm::max(max, point);
    }

    inline bool Hit(const Ray& ray, float min_t, float max_t) const{// add interval
        return HitDistance(ray, min_t, max_t) < std::numeric_limits<float>::infinity();
    }

    inline bool Hit(const Ray& ray, float max_t) const{// add interval
        return HitDistance(ray, max_t) < std::numeric_limits<float>::infinity();
    }

    inline float HitDistance(const Ray& ray, float min_t, float max_t) const{// add interval
#if defined(__SSE__) || defined(_M_AMD64) || defined(_M_X64)
        const __m128 bmin4 = _mm_set_ps(0, min.z, min.y, min.x);
        const __m128 bmax4 = _mm_set_ps(0, max.z, max.y, max.x);

        const __m128 t1 = _mm_mul_ps(_mm_sub_ps(bmin4, ray.O4), ray.rD4);
        const __m128 t2 = _mm_mul_ps(_mm_sub_ps(bmax4, ray.O4), ray.rD4);
        const __m128 vmax4 = _mm_max_ps(t1, t2);
        const __m128 vmin4 = _mm_min_ps(t1, t2);

        alignas(16) float temp_min[4], temp_max[4];
        _mm_store_ps(temp_min, vmin4);
        _mm_store_ps(temp_max, vmax4);
        //__m128 tmp = _mm_max_ps(vmin4, _mm_movehl_ps(vmin4,vmin4));// max ( [0,z,0,z] [0,z,y,x] ) -> max of z and x
        //tmp = _mm_max_ss(tmp, _mm_shuffle_ps(tmp,tmp,_MM_SHUFFLE(1,0,3,2)));// max ( maxzx )?
        //const float tEntry = _mm_cvtss_f32(tmp);
        const float tEntry = std::max(std::max(temp_min[0], temp_min[1]), temp_min[2]);
        const float tExit = std::min(std::min(temp_max[0], temp_max[1]), temp_max[2]);
#else       
        glm::vec3 invD = ray.inv_dir;//store in ray!!!! test speed
        // Calculate t for each slab
        glm::vec3 t0s = (min - ray.origin) * invD;
        glm::vec3 t1s = (max - ray.origin) * invD;
        // Component‐wise min/max
        glm::vec3 tSmaller = glm::min(t0s, t1s);
        glm::vec3 tLarger = glm::max(t0s, t1s);
        // Over all axes
        float tEntry = std::max(std::max(tSmaller.x, tSmaller.y), tSmaller.z);
        float tExit = std::min(std::min(tLarger.x, tLarger.y), tLarger.z);
#endif
        return (tEntry <= tExit && tEntry <= max_t && tExit >= min_t) ? tEntry : std::numeric_limits<float>::infinity();
    }

    inline float HitDistance(const Ray& ray, float max_t) const{// add interval
#if defined(__SSE__) || defined(_M_AMD64) || defined(_M_X64)
        const __m128 bmin4 = _mm_set_ps(0, min.z, min.y, min.x);
        const __m128 bmax4 = _mm_set_ps(0, max.z, max.y, max.x);

        const __m128 t1 = _mm_mul_ps(_mm_sub_ps(bmin4, ray.O4), ray.rD4);
        const __m128 t2 = _mm_mul_ps(_mm_sub_ps(bmax4, ray.O4), ray.rD4);
        __m128 vmax4 = _mm_max_ps(t1, t2);
        __m128 vmin4 = _mm_min_ps(t1, t2);

        alignas(16) float temp_min[4], temp_max[4];
        _mm_store_ps(temp_min, vmin4);
        _mm_store_ps(temp_max, vmax4);
        /*
        __m128 a = _mm_shuffle_ps(vmin4,vmin4,_MM_SHUFFLE(0,3,2,1));//y
        __m128 b = _mm_shuffle_ps(vmin4,vmin4,_MM_SHUFFLE(1,0,3,2));//z
        vmin4 = _mm_max_ps(vmin4,_mm_max_ps(a,b));
        const float tEntry = _mm_cvtss_f32(vmin4);

        a = _mm_shuffle_ps(vmax4,vmax4,_MM_SHUFFLE(0,3,2,1));//y
        b = _mm_shuffle_ps(vmax4,vmax4,_MM_SHUFFLE(1,0,3,2));//z
        vmax4 = _mm_min_ps(vmax4,_mm_min_ps(a,b));
        const float tExit = _mm_cvtss_f32(vmax4);
        */
        const float tEntry = std::max(std::max(temp_min[0], temp_min[1]), temp_min[2]);
        const float tExit = std::min(std::min(temp_max[0], temp_max[1]), temp_max[2]);
#else       
        glm::vec3 invD = ray.inv_dir;//store in ray!!!! test speed
        // Calculate t for each slab
        glm::vec3 t0s = (min - ray.origin) * invD;
        glm::vec3 t1s = (max - ray.origin) * invD;
        // Component‐wise min/max
        glm::vec3 tSmaller = glm::min(t0s, t1s);
        glm::vec3 tLarger = glm::max(t0s, t1s);
        // Over all axes
        float tEntry = std::max(std::max(tSmaller.x, tSmaller.y), tSmaller.z);
        float tExit = std::min(std::min(tLarger.x, tLarger.y), tLarger.z);
#endif
        return (tEntry <= tExit && tExit >= shadowEpsilon && tEntry <= max_t) ? tEntry : std::numeric_limits<float>::infinity();//maybe texit greater than 0 ?
    }


    //using AABB::nothing
    //AABB::universe
    //mediumintr takes phasefunc

    inline float Area() const{
        glm::vec3 e = max - min; // box extent
        return e.x * e.y + e.y * e.z + e.z * e.x;
    }

    glm::vec3 operator[](int i) const{ return (i == 0) ? min : max; }
    glm::vec3& operator[](int i){ return (i == 0) ? min : max; }

    glm::vec3 Corner(int n){
        return { (*this)[(n & 1)].x,
                (*this)[(n & 2) ? 1 : 0].y,
                (*this)[(n & 4) ? 1 : 0].z };
    }
};

struct Bounds2i{
    glm::ivec2 min = { 0,0 };
    glm::ivec2 max = { 0,0 };
};








