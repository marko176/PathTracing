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
    return light->Power() / totalPower;
}
glm::vec3 PowerLightSampler::SampleLd(const Ray& curr_ray,const SurfaceInteraction& interaction,const TLAS& bvh,float u,glm::vec2& UV) const {
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
void PowerLightSampler::PreProcess(const AABB& bbox){
    lightPowers.reserve(lights.size());
    for(const std::shared_ptr<Light>& light : lights){
        light->PreProcess(bbox);
        lightPowers.push_back(light->Power());
        totalPower += lightPowers.back();
    }
}


void RISLightSampler::Add(const std::shared_ptr<Light>& light) {
    lights.push_back(light);
}
void RISLightSampler::Add(const std::vector<std::shared_ptr<Light>>& lights) {
    this->lights.insert(this->lights.end(),lights.begin(),lights.end());
}
   
std::shared_ptr<Light> RISLightSampler::Sample(float u) const {
    //no op
    return nullptr;
}
float RISLightSampler::PMF(const std::shared_ptr<Light>& light) const {
    //no op
    return std::numeric_limits<float>::infinity();
}
glm::vec3 RISLightSampler::SampleLd(const Ray& curr_ray,const SurfaceInteraction& interaction,const TLAS& bvh,float u,glm::vec2& UV) const {
    std::vector<std::pair<glm::vec2,double>> powers;
    powers.reserve(lights.size()*samplesPerLight);
    double totalWeight = 0;
    //test out with weight = 1 for everything -> should be correct
    for(auto& light : lights){
        
        for(int i = 0;i<samplesPerLight;i++){
            glm::vec2 uv = {random_double(),random_double()};
            LightSample lightSample = light->sample(uv);
            glm::vec3 to_light = lightSample.interaction.p - interaction.p;
            Ray shadow_ray(interaction.p, glm::normalize(to_light));
            glm::vec3 L = light->L(lightSample.interaction,shadow_ray);//fix this {}
            glm::vec3 D = interaction.mat->calc_attenuation(curr_ray,interaction,shadow_ray);//has cos term
            powers.push_back({uv,glm::dot(interaction.ns,shadow_ray.dir) > 0 ? glm::dot(L * D,glm::vec3(0.2126f, 0.7152f, 0.0722f)) / light->PDF(lightSample.interaction,shadow_ray): 0});
            totalWeight += powers.back().second;
        }
    }
    if(totalWeight<=0)return {0,0,0};
    glm::vec3 acc = {0,0,0};
    for(int k = 0;k<N;k++){
        float currPower = 0;
        float power = random_double()*totalWeight;
        for(int i = 0;i<powers.size();i++){
            currPower+=powers[i].second;
            if(currPower >= power){
                LightSample lightSample = lights[i/samplesPerLight]->sample(powers[i].first);
                glm::vec3 to_light = lightSample.interaction.p - interaction.p;
                Ray shadow_ray(interaction.p, glm::normalize(to_light));
                if(!bvh.IntersectPred(shadow_ray,glm::length(to_light)-0.005f)){
                    float light_pdf = (powers[i].second / totalWeight) * samplesPerLight * lights[i/samplesPerLight]->PDF(lightSample.interaction,shadow_ray);
                    if(light_pdf <= 0)return {0,0,0};
                    float w2 = light_pdf*light_pdf;
                    float w1 = interaction.mat->PDF(curr_ray,interaction,shadow_ray);
                    w1 = w1*w1;
                    float w_light = (w2) / (w1 + w2);
                    acc+= lights[i/samplesPerLight]->L(lightSample.interaction,shadow_ray) * interaction.mat->calc_attenuation(curr_ray,interaction,shadow_ray) * w_light / light_pdf;
                }
                break;
            }
        }

    }
    return acc/float(N);
}