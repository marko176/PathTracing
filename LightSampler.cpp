#include "LightSampler.hpp"

void UniformLightSampler::Add(const std::shared_ptr<Light>& light) {
    lights.push_back(light);
}
void UniformLightSampler::Add(const std::vector<std::shared_ptr<Light>>& lights) {
    this->lights.insert(this->lights.end(),lights.begin(),lights.end());
}
std::shared_ptr<Light> UniformLightSampler::Sample(float u) const {
    int index = std::min<int>((u * lights.size()),lights.size()-1);
    return lights[index];
}
float UniformLightSampler::PMF(const std::shared_ptr<Light>& light) const {
    return 1.0f / lights.size();
}
glm::vec3 UniformLightSampler::SampleLd(const Ray& curr_ray,const SurfaceInteraction& interaction,const TLAS& bvh,float u,glm::vec2& UV) const {
    std::shared_ptr<Light> sampled_light = Sample(u);
    if(sampled_light == nullptr)return {0,0,0};
    if(sampled_light->isDelta()){
        glm::vec3 dir = sampled_light->sample(UV).dir;
        Ray shadow_ray(interaction.p,dir);
        if(glm::dot(interaction.ns,shadow_ray.dir) > 0 && !bvh.IntersectPred(shadow_ray,1e30f)){
            return sampled_light->L({},shadow_ray) * interaction.mat->calc_attenuation(curr_ray,interaction,shadow_ray) / PMF(sampled_light);
        }
    }else{
        LightSample lightSample = sampled_light->sample(UV);
        glm::vec3 to_light = lightSample.interaction.p - interaction.p;
        Ray shadow_ray(interaction.p, glm::normalize(to_light));
        if(glm::dot(interaction.ns,shadow_ray.dir) > 0 && !bvh.IntersectPred(shadow_ray,glm::length(to_light)-0.005f)){
            float light_pdf = PMF(sampled_light) * sampled_light->PDF(lightSample.interaction,shadow_ray);
            if(light_pdf <= 0)return {0,0,0};
            float w2 = light_pdf*light_pdf;
            float w1 = interaction.mat->PDF(curr_ray,interaction,shadow_ray);
            w1 = w1*w1;
            float w_light = (w2) / (w1 + w2);
            return sampled_light->L(lightSample.interaction,shadow_ray) * interaction.mat->calc_attenuation(curr_ray,interaction,shadow_ray) * w_light / light_pdf;
        }
    }
    return {0,0,0};
}


void PowerLightSampler::Add(const std::shared_ptr<Light>& light) {
    lights.push_back(light);
}
void PowerLightSampler::Add(const std::vector<std::shared_ptr<Light>>& lights) {
    this->lights.insert(this->lights.end(),lights.begin(),lights.end());
}
   
std::shared_ptr<Light> PowerLightSampler::Sample(float u) const {
    if(lights.empty())return nullptr;
    float currPower = 0;
    float power = u*totalPower;
    for(int i = 0;i<lightPowers.size();i++){
        currPower+=lightPowers[i];
        if(currPower >= power){
            return lights[i];
        }
    }
    return lights.back();
}
float PowerLightSampler::PMF(const std::shared_ptr<Light>& light) const {
    if(totalPower==0)return 1.0f;
    return light->Power() / totalPower;
}

inline bool IntersectTr(Ray ray, float t, const TLAS& bvh, glm::vec3& Tr){
    Tr = {1,1,1};
    while(true){
        SurfaceInteraction interaction;
        bool hit = bvh.Intersect(ray,interaction,t);
        if(ray.medium)
            Tr *= ray.medium->Tr(ray,interaction.t);
        
        if(!hit)
            return false;
        if(interaction.mat != nullptr)
            return true;

        ray = Ray(ray.at(interaction.t),ray.dir,interaction.getMedium(ray.dir));
        t-=interaction.t;
        
    }
    return false;
}

glm::vec3 PowerLightSampler::SampleLd(const Ray& curr_ray,const SurfaceInteraction& interaction,const TLAS& bvh,float u,glm::vec2& UV) const {
    std::shared_ptr<Light> sampled_light = Sample(u);
    if(sampled_light == nullptr)return {0,0,0};
    glm::vec3 Tr = {1,1,1};
    if(sampled_light->isDelta()){
        glm::vec3 dir = sampled_light->sample(UV).dir;
        Ray shadow_ray(interaction.p,dir);
        if(glm::dot(interaction.ns,shadow_ray.dir) > 0 && !bvh.IntersectPred(shadow_ray,1e30f)){
            float light_pdf = PMF(sampled_light);
            if(light_pdf <= 0)return {0,0,0};
            return Tr * sampled_light->L({},shadow_ray) * interaction.mat->calc_attenuation(curr_ray,interaction,shadow_ray) / light_pdf;
        }
    }else{
        LightSample lightSample = sampled_light->sample(UV);
        glm::vec3 to_light = lightSample.interaction.p - interaction.p;
        Ray shadow_ray(interaction.p, glm::normalize(to_light));
        if(glm::dot(interaction.ns,shadow_ray.dir) > 0 && !bvh.IntersectPred(shadow_ray,glm::length(to_light)-0.005f)/*!bvh.IntersectPred(shadow_ray,glm::length(to_light)-0.005f)*/){
            float light_pdf = PMF(sampled_light) * sampled_light->PDF(lightSample.interaction,shadow_ray);
            if(light_pdf <= 0)return {0,0,0};
            float w2 = light_pdf*light_pdf;
            float w1 = interaction.mat->PDF(curr_ray,interaction,shadow_ray);
            w1 = w1*w1;
            float w_light = (w2) / (w1 + w2);
            return Tr * sampled_light->L(lightSample.interaction,shadow_ray) * interaction.mat->calc_attenuation(curr_ray,interaction,shadow_ray) * w_light / light_pdf;
        }
 
    }
    return {0,0,0};
}
void PowerLightSampler::PreProcess(const AABB& bbox){
    lightPowers.reserve(lights.size());
    for(const std::shared_ptr<Light>& light : lights){
        light->PreProcess(bbox);
        lightPowers.push_back(light->Power());
        totalPower += lightPowers.back();
    }
}


