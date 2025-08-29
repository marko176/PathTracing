#include "LightSampler.hpp"

void UniformLightSampler::Add(Light* light) {
    lights.push_back(light);
}
void UniformLightSampler::Add(const std::vector<Light*>& lights) {
    this->lights.insert(this->lights.end(),lights.begin(),lights.end());
}
const Light* UniformLightSampler::Sample(float u) const {
    int index = std::min<int>((u * lights.size()),lights.size()-1);
    return lights[index];
}
float UniformLightSampler::PMF(const Light* light) const {
    return 1.0f / lights.size();
}
glm::vec3 UniformLightSampler::SampleLd(const Ray& curr_ray,const SurfaceInteraction& interaction,float brdf,const TLAS& bvh,float u,glm::vec2& UV) const {
    const Light* sampled_light = Sample(u);
    if(sampled_light->isDelta()){
        glm::vec3 dir = sampled_light->sample(UV).dir;
        Ray shadow_ray(interaction.p,dir);
        if(glm::dot(interaction.n,shadow_ray.dir) > 0 && !bvh.IntersectPred(shadow_ray,1e30f)){
            return sampled_light->L({},shadow_ray) * interaction.mat->calc_attenuation(curr_ray,interaction,shadow_ray) / PMF(sampled_light);
        }
    }else{
        LightSample lightSample = sampled_light->sample(UV);
        glm::vec3 to_light = lightSample.interaction.p - interaction.p;
        Ray shadow_ray(interaction.p, glm::normalize(to_light));
        if(glm::dot(interaction.n,shadow_ray.dir) > 0 && !bvh.IntersectPred(shadow_ray,glm::length(to_light)-0.005f)){
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


void PowerLightSampler::Add(Light* light) {
    lights.push_back(light);
}
void PowerLightSampler::Add(const std::vector<Light*>& lights) {
    this->lights.insert(this->lights.end(),lights.begin(),lights.end());
}
   
const Light* PowerLightSampler::Sample(float u) const {
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
float PowerLightSampler::PMF(const Light* light) const {
    return light->Power() / totalPower;
}
glm::vec3 PowerLightSampler::SampleLd(const Ray& curr_ray,const SurfaceInteraction& interaction,float brdf,const TLAS& bvh,float u,glm::vec2& UV) const {
    const Light* sampled_light = Sample(u);
    if(sampled_light->isDelta()){
        glm::vec3 dir = sampled_light->sample(UV).dir;
        Ray shadow_ray(interaction.p,dir);
        if(glm::dot(interaction.n,shadow_ray.dir) > 0 && !bvh.IntersectPred(shadow_ray,1e30f)){
            return sampled_light->L({},shadow_ray) * interaction.mat->calc_attenuation(curr_ray,interaction,shadow_ray) / PMF(sampled_light);
        }
    }else{
        LightSample lightSample = sampled_light->sample(UV);
        glm::vec3 to_light = lightSample.interaction.p - interaction.p;
        Ray shadow_ray(interaction.p, glm::normalize(to_light));
        if(glm::dot(interaction.n,shadow_ray.dir) > 0 && !bvh.IntersectPred(shadow_ray,glm::length(to_light)-0.005f)){
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
void PowerLightSampler::PreProcess(const AABB& bbox){
    lightPowers.reserve(lights.size());
    for(Light* light : lights){
        light->PreProcess(bbox);
        lightPowers.push_back(light->Power());
        totalPower += lightPowers.back();
    }
}