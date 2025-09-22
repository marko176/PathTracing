#pragma once
#include "Film.hpp"
class Camera {
public: 
    Camera(const glm::dvec3& lookFrom, const glm::dvec3& lookAt, double fov,const std::shared_ptr<Film>& film) : Camera(lookFrom,lookAt,fov,film,0,0) {

    };
    Camera(const glm::dvec3& lookFrom, const glm::dvec3& lookAt, double fov,const std::shared_ptr<Film>& film, double FocusAngle,double FocusDistance) : lookFrom(lookFrom), lookAt(lookAt) , Fov(fov) , film(film), FocusAngle(FocusAngle), FocusDistance(FocusDistance) {
        w = glm::normalize(lookFrom-lookAt);
        u = glm::normalize(glm::cross({0,1,0},w));//up is {0,1,0}
        v = glm::cross(w,u);
        defocusRadius = FocusDistance * std::tan(FocusAngle/2.0);
        halfWidth  = std::tan(fov * 0.5);
        halfHeight = halfWidth * film->Resolution().y / film->Resolution().x;
    };

    Ray GenerateRay(const glm::vec2& p,const glm::vec2& LensUV = {0,0}, float time = 0) const {
        double u_coord = p.x / double(film->Resolution().x);
        double v_coord = p.y / double(film->Resolution().y);
        glm::dvec3 direction = glm::normalize(-w + (2.0f * u_coord - 1.0f) * halfWidth * u + (2.0f * v_coord - 1.0f) * halfHeight * v);
        if(FocusDistance == 0 || FocusAngle == 0){
            return Ray(lookFrom,direction);
        }
        glm::vec2 pLens = inUnitDisk(LensUV);
        glm::dvec3 defocus_disk_u = u * defocusRadius;
        glm::dvec3 defocus_disk_v = v * defocusRadius;
        direction = direction*FocusDistance;
        glm::dvec3 offset = glm::dvec3(pLens.x) * defocus_disk_u + glm::dvec3(pLens.y) * defocus_disk_v;
        return Ray(lookFrom + offset,glm::normalize(direction - offset));
    }

    std::shared_ptr<Film> GetFilm() const {
        return film;
    }
protected:
    glm::dvec3 lookFrom;
    glm::dvec3 lookAt;
    double Fov;
    std::shared_ptr<Film> film;
    double FocusAngle;
    double FocusDistance;
    glm::dvec3 w;
    glm::dvec3 u;
    glm::dvec3 v;
    double defocusRadius;
    double halfWidth;
    double halfHeight;
};