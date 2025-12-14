#pragma once 
#include "Interaction.hpp"
class onb{
public:
    onb(const glm::vec3& n){
        axis[2] = n;
        glm::vec3 up = (std::fabs(axis[2].x) > 0.9999) ? glm::vec3(0, 1, 0) : glm::vec3(1, 0, 0);
        axis[1] = glm::normalize(glm::cross(axis[2], up));
        axis[0] = glm::cross(axis[1], axis[2]);
    }

    onb(const SurfaceInteraction& interaction){
        axis[2] = interaction.ns;
        axis[0] = interaction.tangent;
        axis[1] = glm::cross(axis[2], axis[0]);
    }

    constexpr glm::vec3 toWorld(const glm::vec3& v) const{
        return v.x * axis[0] + v.y * axis[1] + v.z * axis[2];
    }

    constexpr glm::vec3 toLocal(const glm::vec3& v) const{
        return glm::vec3(glm::dot(v, axis[0]),
                        glm::dot(v, axis[1]),
                        glm::dot(v, axis[2]));
    }

private:
    glm::vec3 axis[3];
};