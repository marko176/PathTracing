#pragma once
#include "Light.hpp"
#include "Primitive.hpp"
class LightSampler {
public:
    virtual void Add(Light* light) = 0;
    virtual void Add(const std::vector<Light*>& lights) {
        for(Light* l : lights){
            Add(l);
        }
    }
    virtual const Light* Sample(float u) const = 0;
    virtual float PMF(const Light* light) const = 0;
    virtual glm::vec3 SampleLd(const Ray& curr_ray,const SurfaceInteraction& interaction,float brdf,const TLAS& bvh,float u,glm::vec2& UV) const = 0;
    virtual void PreProcess(const AABB& bbox) {}
};

class UniformLightSampler : public LightSampler {
public:
    void Add(Light* light) override ;

    virtual void Add(const std::vector<Light*>& lights) ;

    const Light* Sample(float u) const override ;

    float PMF(const Light* light) const override ;

    glm::vec3 SampleLd(const Ray& curr_ray,const SurfaceInteraction& interaction,float brdf,const TLAS& bvh,float u,glm::vec2& UV) const override ;
private:
    std::vector<Light*> lights;
};

class PowerLightSampler : public LightSampler {
public: 
    void Add(Light* light) override ;

    virtual void Add(const std::vector<Light*>& lights) ;
       
    const Light* Sample(float u) const override ;

    float PMF(const Light* light) const override ;

    glm::vec3 SampleLd(const Ray& curr_ray,const SurfaceInteraction& interaction,float brdf,const TLAS& bvh,float u,glm::vec2& UV) const override ;

    void PreProcess(const AABB& bbox) override ;
private:
    std::vector<Light*> lights;
    std::vector<float> lightPowers;
    float totalPower;
};
