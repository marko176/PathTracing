#pragma once
#include "Interaction.hpp"
#include "LightSampler.hpp"
#include "Util.hpp"
class Scene;
class Camera;
class Sampler;


class Integrator {
public:
    virtual ~Integrator() = default;
    virtual void Render() const = 0;
    Integrator(const std::shared_ptr<Scene>& scene, const std::shared_ptr<Camera>& camera,const std::shared_ptr<Sampler>& sampler) : scene(scene), camera(camera), sampler(sampler) {}
    bool Unoccluded(const Ray& ray, float t) const;
    virtual glm::vec3 Li(Ray ray) const = 0;
    bool Intersect(const Ray& ray, SurfaceInteraction& interaction, float max = std::numeric_limits<float>::infinity()) const ;
    bool IntersectTr(const Ray& ray, SurfaceInteraction& interaction, glm::vec3& Tr, float max = std::numeric_limits<float>::infinity()) const;
protected:
    std::shared_ptr<Scene> scene;
    std::shared_ptr<Camera> camera;
    std::shared_ptr<Sampler> sampler;
};

class TileIntegrator : public Integrator {
public:
    virtual ~TileIntegrator() = default;
    TileIntegrator(const std::shared_ptr<Scene>& scene, const std::shared_ptr<Camera>& camera,const std::shared_ptr<Sampler>& sampler) : Integrator(scene,camera,sampler) {}
    
    void Render() const override ;
};

class PathIntegrator : public TileIntegrator {
public:
    virtual ~PathIntegrator() = default;
    PathIntegrator(const std::shared_ptr<Scene>& scene, const std::shared_ptr<Camera>& camera,const std::shared_ptr<Sampler>& sampler, int maxDepth) : PathIntegrator(scene,camera,sampler,std::make_shared<UniformLightSampler>(),maxDepth) {}
    PathIntegrator(const std::shared_ptr<Scene>& scene, const std::shared_ptr<Camera>& camera,const std::shared_ptr<Sampler>& sampler, const std::shared_ptr<LightSampler>& lightSampler, int maxDepth) : TileIntegrator(scene,camera,sampler), lightSampler(lightSampler), maxDepth(maxDepth) {}
    
    glm::vec3 Li(Ray ray) const override ;
protected:
    glm::vec3 SampleLd(const Ray& ray,const SurfaceInteraction& interaction,float u,const glm::vec2& UV) const ;

    std::shared_ptr<LightSampler> lightSampler;
    uint32_t maxDepth;
};

class VolPathIntegrator : public TileIntegrator {
public:
    virtual ~VolPathIntegrator() = default;
    VolPathIntegrator(const std::shared_ptr<Scene>& scene, const std::shared_ptr<Camera>& camera,const std::shared_ptr<Sampler>& sampler, int maxDepth) : VolPathIntegrator(scene,camera,sampler,std::make_shared<UniformLightSampler>(),maxDepth) {}
    VolPathIntegrator(const std::shared_ptr<Scene>& scene, const std::shared_ptr<Camera>& camera,const std::shared_ptr<Sampler>& sampler, const std::shared_ptr<LightSampler>& lightSampler, int maxDepth) : TileIntegrator(scene,camera,sampler), lightSampler(lightSampler), maxDepth(maxDepth) {}
    
    glm::vec3 Li(Ray ray) const override ;
protected:
    glm::vec3 SampleLd(const Ray& ray,const SurfaceInteraction& interaction,float u,const glm::vec2& UV) const ;
    glm::vec3 SampleLdMedium(const Ray& ray,const MediumInteraction& interaction,float u,const glm::vec2& UV) const ;
    
    std::shared_ptr<LightSampler> lightSampler;
    uint32_t maxDepth;
};
