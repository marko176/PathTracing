#pragma once
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtx/intersect.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <memory>
class Medium;
struct Ray{
    glm::vec3 origin = {0,0,0};
    glm::vec3 dir = {0,0,0};
    glm::vec3 inv_dir = {0,0,0};
    std::shared_ptr<Medium> medium = nullptr;
    float time = 0;
    Ray() = default;
    Ray(const glm::vec3& origin, const glm::vec3& dir,const std::shared_ptr<Medium>& medium = nullptr) : origin(origin), dir(dir), inv_dir(1.0f/dir), medium(medium), time(0) {}
    Ray(const glm::vec3& origin, const glm::vec3& dir,float time, const std::shared_ptr<Medium>& medium = nullptr) : origin(origin), dir(dir), inv_dir(1.0f/dir), medium(medium), time(time) {}

    inline glm::vec3 at(float t) const {
        return origin + t * dir;
    }
};