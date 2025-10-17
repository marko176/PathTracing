#pragma once 
#include "Interaction.hpp"
class onb {
public:
    onb(const glm::vec3& n) {
        axis[2] = glm::normalize(n);
        glm::vec3 a = (std::fabs(axis[2].x) > 0.9999) ? glm::vec3(0,1,0) : glm::vec3(1,0,0);
        axis[1] = glm::normalize(glm::cross(axis[2], a));
        axis[0] = glm::cross(axis[1],axis[2]);
    }

    onb(const SurfaceInteraction& interaction) {
        axis[2] = glm::normalize(interaction.ns);
        if(interaction.tangent == glm::vec3(0,0,0)){
            glm::vec3 a = (std::fabs(axis[2].x) > 0.9999) ? glm::vec3(0,1,0) : glm::vec3(1,0,0);
            axis[1] = glm::normalize(glm::cross(axis[2], a));
            axis[0] = glm::cross(axis[2],axis[1]);
        }else{    
            axis[0] = glm::normalize(interaction.tangent - axis[2] * glm::dot(axis[2],interaction.tangent));     // T 
            axis[1] = glm::cross(axis[2],axis[0]);
            if(glm::dot(axis[1],interaction.bitangent)<0){
                axis[0]=-axis[0];
                axis[1]=-axis[1];
            }

        }
    }

    const glm::vec3& u() const { return axis[0]; }
    const glm::vec3& v() const { return axis[1]; }
    const glm::vec3& w() const { return axis[2]; }

    glm::vec3 toWorld(const glm::vec3& v) const {
        // v = (v.x, v.y, v.z) in tangent space => world = v.x * T + v.y * B + v.z * N
        return v.x * axis[0] + v.y * axis[1] + v.z * axis[2];
    }

    // world→local
    glm::vec3 toLocal(const glm::vec3& v) const {
        // inverse of orthonormal matrix is its transpose:
        return glm::vec3(glm::dot(v, axis[0]),
                         glm::dot(v, axis[1]),
                         glm::dot(v, axis[2]));
    }

private:
    glm::vec3 axis[3];
};