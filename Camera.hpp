#pragma once
#include "Film.hpp"
#include "Medium.hpp"
class Camera{
public:
    Camera(const glm::vec3& lookFrom, const glm::vec3& lookAt, float fov, const std::shared_ptr<Film>& film) : Camera(lookFrom, lookAt, fov, film, 0, 0){};
    Camera(const glm::vec3& lookFrom, const glm::vec3& lookAt, float fov, const std::shared_ptr<Film>& film, float FocusAngle, float FocusDistance, const std::shared_ptr<Medium>& medium = nullptr) : lookFrom(lookFrom), lookAt(lookAt), Fov(fov), film(film), FocusAngle(FocusAngle), FocusDistance(FocusDistance), cameraMedium(medium){
        w = glm::normalize(lookFrom - lookAt);
        u = glm::normalize(glm::cross({ 0,1,0 }, w));//up is {0,1,0}
        v = glm::cross(w, u);
        defocusRadius = FocusDistance * std::tan(FocusAngle / 2.0);
        halfWidth = std::tan(fov * 0.5);
        halfHeight = halfWidth * film->Resolution().y / film->Resolution().x;
    };

    Camera(const glm::vec3& lookFrom, const glm::vec3& lookAt, float fov, const std::shared_ptr<Film>& film, const glm::vec2& shutterBounds) : Camera(lookFrom, lookAt, fov, film, 0, 0){
        shutterStart = shutterBounds.x;
        shutterEnd = shutterBounds.y;
    };

    Ray GenerateRay(const glm::vec2& p, float time, const glm::vec2& LensUV = { 0,0 }) const{
        float u_coord = p.x / static_cast<float>(film->Resolution().x);
        float v_coord = p.y / static_cast<float>(film->Resolution().y);
        glm::vec3 direction = glm::normalize(-w + (2.0f * u_coord - 1.0f) * halfWidth * u + (2.0f * v_coord - 1.0f) * halfHeight * v);
        float t = glm::mix(shutterStart, shutterEnd, time);
        if(FocusDistance == 0 || FocusAngle == 0){
            return Ray(lookFrom, direction, t, GetMedium());
        }
        glm::vec2 pLens = inUnitDisk(LensUV);
        glm::vec3 defocus_disk_u = u * defocusRadius;
        glm::vec3 defocus_disk_v = v * defocusRadius;
        direction *= FocusDistance;
        glm::vec3 offset = glm::vec3(pLens.x) * defocus_disk_u + glm::vec3(pLens.y) * defocus_disk_v;
        return Ray(lookFrom + offset, glm::normalize(direction - offset), t, GetMedium());
    }

    std::shared_ptr<Film> GetFilm() const{
        return film;
    }

    std::shared_ptr<Medium> GetMedium() const{
        return cameraMedium;
    }

    void SetMedium(const std::shared_ptr<Medium>& medium){
        cameraMedium = medium;
    }
protected:
    glm::vec3 lookFrom;
    glm::vec3 lookAt;
    float Fov;
    std::shared_ptr<Film> film;
    float FocusAngle;
    float FocusDistance;
    glm::vec3 w;
    glm::vec3 u;
    glm::vec3 v;
    float defocusRadius;
    float halfWidth;
    float halfHeight;
    std::shared_ptr<Medium> cameraMedium;
    float shutterStart;
    float shutterEnd;
};