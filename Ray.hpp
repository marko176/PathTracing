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
    Ray(const glm::vec3& origin, const glm::vec3& dir,const std::shared_ptr<Medium>& medium = nullptr) : origin(origin), inv_dir(1.0f/dir), dir(dir), medium(medium), time(0) {

    }
    Ray(const glm::vec3& origin, const glm::vec3& dir,float time, const std::shared_ptr<Medium>& medium = nullptr) : origin(origin), inv_dir(1.0f/dir),dir(dir), medium(medium), time(time) {

    }

    inline glm::vec3 at(float t) const {
        return origin + t * dir;
    }
};