#pragma once
#include "Ray.hpp"
#include "Interaction.hpp"
#include "PhaseFunction.hpp"
class Medium{
public:
    virtual glm::vec3 Tr(const Ray& ray, float t) const = 0;
    virtual glm::vec3 Sample(const Ray& ray, float t, MediumInteraction& interaction) const = 0;
    virtual bool IsEmmisive() const = 0;
    virtual glm::vec3 Le() const = 0;
    virtual ~Medium() = default;
};

class HomogeneusMedium : public Medium{
public:
    virtual ~HomogeneusMedium() = default;

    HomogeneusMedium(const glm::vec3& sigma_a, const glm::vec3& sigma_s, const std::shared_ptr<PhaseFunction> phaseFunction, float density = 1.0f, const glm::vec3& Le = glm::vec3 { 0,0,0 }, float LeDensity = 1.0f) : sigma_a(density* sigma_a), sigma_s(density* sigma_s), sigma_t(density* (sigma_a + sigma_s)), emmision(Le* LeDensity), phaseFunction(phaseFunction){

    }

    glm::vec3 Tr(const Ray& ray, float t) const override{
        //t should be ray.t!!
        return glm::exp(-sigma_t * std::min(t, std::numeric_limits<float>::max()));
    }

    glm::vec3 Sample(const Ray& ray, float t, MediumInteraction& interaction) const override{
        int channel = (int)random_float(0, 3);

        float scatterDist = std::min<float>(-std::log(1.0 - random_float()) / sigma_t[channel], t);//use sampeler variable
        bool sampledMedium = scatterDist < t;
        if(sampledMedium){
            //enable_shared_from this
            interaction = MediumInteraction(ray.at(scatterDist), { 0,0,0 }, ray.medium, phaseFunction);//should set medium to itself but it is shared ptr?
        }

        glm::vec3 tr = Tr(ray, scatterDist);

        glm::vec3 density = sampledMedium ? (sigma_t * tr) : tr;
        float pdf = 0;
        for(int i = 0;i < 3;i++){
            pdf += density[i];
        }
        pdf /= 3.0;
        return sampledMedium ? (tr * sigma_s / pdf) : (tr / pdf);
    }

    bool IsEmmisive() const override{
        return Le() != glm::vec3(0, 0, 0);
    }

    glm::vec3 Le() const override{
        return emmision;
    }
private:
    glm::vec3 sigma_a;
    glm::vec3 sigma_s;
    glm::vec3 sigma_t;
    glm::vec3 emmision;
    std::shared_ptr<PhaseFunction> phaseFunction;
};