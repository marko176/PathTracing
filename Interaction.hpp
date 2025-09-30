#pragma once
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtx/intersect.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <memory>

class Material;
class Medium;
class PhaseFunction;
class Light;

struct GeometricInteraction {
    glm::vec3 p = {0,0,0};
    glm::vec3 n = {0,0,0};//geometric normal
    std::shared_ptr<Medium> medium = nullptr;
    
    bool isSurfaceInteraction() const {
        return n != glm::vec3(0,0,0);
    }

    bool isMediumInteraction() const {
        return !isSurfaceInteraction();
    }

    std::shared_ptr<Medium> getMedium(const glm::vec3& dir) const{
        if(glm::dot(dir,n)<0)return medium;//going into surface
        return nullptr;
    }
    GeometricInteraction() = default;
    GeometricInteraction(const glm::vec3& point, const glm::vec3& geometric_normal,const std::shared_ptr<Medium>& medium = nullptr) : p{point} , n{geometric_normal} , medium(medium) {

    }
};

struct SurfaceInteraction : public GeometricInteraction{
    SurfaceInteraction(const glm::vec3& point, const glm::vec3& geometric_normal, const glm::vec2& uv,const std::shared_ptr<Medium>& medium = nullptr) : GeometricInteraction{point,geometric_normal,medium}, uv{uv} {

    }

    SurfaceInteraction() = default;
    SurfaceInteraction(const GeometricInteraction& interaction) : GeometricInteraction{interaction} {

    }
    glm::vec2 uv = {0,0};
    float t = 0;//we dont need t if we have point? t = glm::length(ray)
    glm::vec3 ns;
    std::shared_ptr<Material> mat;
    glm::vec3 tangent;
    glm::vec3 bitangent;
    std::shared_ptr<Light> AreaLight;
};

struct MediumInteraction : public GeometricInteraction {
    MediumInteraction(const glm::vec3& point, const glm::vec3& geometric_normal, const std::shared_ptr<Medium>& medium, const std::shared_ptr<PhaseFunction>& phaseFunction) : GeometricInteraction{point,geometric_normal,medium}, phaseFunction{phaseFunction} {

    }

    MediumInteraction() = default;

    bool isValid() const {
        return phaseFunction != nullptr;
    }
    std::shared_ptr<PhaseFunction> phaseFunction;
};