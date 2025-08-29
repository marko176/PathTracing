#pragma once
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtx/intersect.hpp>
#include <glm/gtc/matrix_transform.hpp>



struct Ray{
    glm::vec3 origin;
    glm::vec3 dir;
    glm::vec3 inv_dir;
    Ray() = default;
    Ray(const glm::vec3& origin, const glm::vec3& dir) : origin(origin), dir(dir), inv_dir(1.0f/dir) {}

    inline glm::vec3 at(float t) const {
        return origin + t * dir;
    }
};