#pragma once
#include "Light.hpp"
#include "Primitive.hpp"
class LightSampler {
public:
    virtual ~LightSampler() = default;
    virtual void Add(const std::shared_ptr<Light>& light) = 0;
    virtual void Add(const std::vector<std::shared_ptr<Light>>& lights) {
        for(const std::shared_ptr<Light>& l : lights){
            Add(l);
        }
    }
    virtual std::shared_ptr<Light> Sample(float u) const = 0;
    virtual float PMF(const std::shared_ptr<Light>& light) const = 0;
    virtual glm::vec3 SampleLd(const Ray& curr_ray,const SurfaceInteraction& interaction,const TLAS& bvh,float u,glm::vec2& UV) const = 0;
    virtual void PreProcess(const AABB& bbox) {}
};

class UniformLightSampler : public LightSampler {
public:
    virtual ~UniformLightSampler() = default;
    
    void Add(const std::shared_ptr<Light>& light) override ;

    void Add(const std::vector<std::shared_ptr<Light>>& lights) override;

    std::shared_ptr<Light> Sample(float u) const override ;

    float PMF(const std::shared_ptr<Light>& light) const override ;

    glm::vec3 SampleLd(const Ray& curr_ray,const SurfaceInteraction& interaction,const TLAS& bvh,float u,glm::vec2& UV) const override ;
private:
    std::vector<std::shared_ptr<Light>> lights;
};

class PowerLightSampler : public LightSampler {
public: 
    virtual ~PowerLightSampler() = default;

    void Add(const std::shared_ptr<Light>& light) override ;

    void Add(const std::vector<std::shared_ptr<Light>>& lights) override;
       
    std::shared_ptr<Light> Sample(float u) const override ;

    float PMF(const std::shared_ptr<Light>& light) const override ;

    glm::vec3 SampleLd(const Ray& curr_ray,const SurfaceInteraction& interaction,const TLAS& bvh,float u,glm::vec2& UV) const override ;

    void PreProcess(const AABB& bbox) override ;
private:
    std::vector<std::shared_ptr<Light>> lights;
    std::vector<float> lightPowers;
    float totalPower;
};

