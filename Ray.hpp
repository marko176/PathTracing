#pragma once
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtx/intersect.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <memory>

#include <immintrin.h>

#define SSE_RAY 0


class Medium;
struct Ray{
#if (defined(__SSE__) || defined(_M_AMD64) || defined(_M_X64)) && SSE_RAY
    union{
        glm::vec3 origin;
        alignas(16) __m128 O4;
    };
    union{
        glm::vec3 inv_dir;
        alignas(16) __m128 rD4;
    };
#else
    glm::vec3 origin;
    glm::vec3 inv_dir;
#endif
    glm::vec3 dir = { 0,0,0 };
    std::shared_ptr<Medium> medium = nullptr;
    float time = 0;
    Ray() = default;
    Ray(const glm::vec3& rayOrigin, const glm::vec3& rayDir, const std::shared_ptr<Medium>& rayMedium = nullptr) : origin(rayOrigin), inv_dir { std::abs(rayDir.x) < 1e-32f ? 1e32f : 1.0f / rayDir.x, std::abs(rayDir.y) < 1e-32f ? 1e32f : 1.0f / rayDir.y,std::abs(rayDir.z) < 1e-32f ? 1e32f : 1.0f / rayDir.z }, dir(rayDir), medium(rayMedium), time(0){

    }
    Ray(const glm::vec3& rayOrigin, const glm::vec3& rayDir, float rayTime, const std::shared_ptr<Medium>& rayMedium = nullptr) : origin(rayOrigin), inv_dir { std::abs(rayDir.x) < 1e-32f ? 1e32f : 1.0f / rayDir.x, std::abs(rayDir.y) < 1e-32f ? 1e32f : 1.0f / rayDir.y,std::abs(rayDir.z) < 1e-32f ? 1e32f : 1.0f / rayDir.z }, dir(rayDir), medium(rayMedium), time(rayTime){

    }

    inline glm::vec3 at(float t) const{
        return origin + t * dir;
    }
};