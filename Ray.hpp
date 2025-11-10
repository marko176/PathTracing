#pragma once
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtx/intersect.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <memory>

#include <immintrin.h>

class Medium;
struct Ray{
#if defined(__SSE__)
    union {
        glm::vec3 origin;
        alignas(16) __m128 O4;
    };
    union {
        glm::vec3 inv_dir;
        alignas(16) __m128 rD4;
    };
#elif
    glm::vec3 origin;
    glm::vec3 inv_dir;
#endif
    glm::vec3 dir = {0,0,0};
    std::shared_ptr<Medium> medium = nullptr;
    float time = 0;
    Ray() = default;
    Ray(const glm::vec3& rayOrigin, const glm::vec3& rayDir,const std::shared_ptr<Medium>& rayMedium = nullptr) : origin(rayOrigin), inv_dir{rayDir.x == 0 ? 1e32f : 1.0f/rayDir.x, rayDir.y == 0 ? 1e32f : 1.0f/rayDir.y,rayDir.z == 0 ? 1e32f : 1.0f/rayDir.z}, dir(rayDir), medium(rayMedium), time(0) {

    }
    Ray(const glm::vec3& rayOrigin, const glm::vec3& rayDir,float rayTime, const std::shared_ptr<Medium>& rayMedium = nullptr) : origin(rayOrigin), inv_dir{rayDir.x == 0 ? 1e32f : 1.0f/rayDir.x, rayDir.y == 0 ? 1e32f : 1.0f/rayDir.y,rayDir.z == 0 ? 1e32f : 1.0f/rayDir.z},dir(rayDir), medium(rayMedium), time(rayTime) {

    }

    inline glm::vec3 at(float t) const {
        return origin + t * dir;
    }
};