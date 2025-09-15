#pragma once
#include "Ray.hpp"
#include "Hit_record.hpp"
class Medium{
public:
    virtual glm::vec3 Tr(const Ray& ray, float t) const = 0;
    virtual glm::vec3 Scatter(const Ray& ray, float t, MediumInteraction& interaction) const = 0;
    virtual ~Medium() = default;
};

class HomogeneusMedium : public Medium {
public:
    HomogeneusMedium(const glm::vec3& sigma_a, const glm::vec3& sigma_s,float density) : sigma_a(sigma_a) , sigma_s(sigma_s), sigma_tr(sigma_a+sigma_s), density(density), negInvDensity(-1.0/density) {

    }

    glm::vec3 Tr(const Ray& ray, float t) const override {
        //t should be ray.t!!
        return glm::exp(density * (-sigma_tr * std::min(t,std::numeric_limits<float>::max())));
    }

    glm::vec3 Scatter(const Ray& ray, float t,MediumInteraction& interaction) const override {


        float scatterDist = negInvDensity * std::log(random_float());//use sampeler variable
        if(t > scatterDist){
            interaction = MediumInteraction(ray.at(scatterDist),{0,0,0},{0,0});
        }else{
            interaction = MediumInteraction(ray.at(scatterDist),{1,1,1},{0,0});
        }
        
        return Tr(ray,t);
    }
private:
    glm::vec3 sigma_a;
    glm::vec3 sigma_s;
    glm::vec3 sigma_tr;
    float density;
    float negInvDensity;
};