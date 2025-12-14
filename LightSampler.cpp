#include "LightSampler.hpp"

void UniformLightSampler::Add(const std::shared_ptr<Light>& light){
    lights.push_back(light);
}

std::shared_ptr<Light> UniformLightSampler::Sample(float u) const{
    if(lights.empty())return nullptr;
    int index = std::min<int>((u * lights.size()), lights.size() - 1);
    return lights[index];
}
float UniformLightSampler::PMF(const std::shared_ptr<Light>& light) const{
    if(lights.empty())return 0;
    return 1.0f / lights.size();
}

void UniformLightSampler::PreProcess(const AABB& bbox){
    std::vector<std::shared_ptr<Light>> validLights;
    validLights.reserve(lights.size());
    for(auto&& light : lights){
        light->PreProcess(bbox);
        if(light->Power() < 0.01f)continue;
        validLights.emplace_back(light);
    }
    lights = validLights;
}


void PowerLightSampler::Add(const std::shared_ptr<Light>& light){
    lights.push_back(light);
}


std::shared_ptr<Light> PowerLightSampler::Sample(float u) const{
    if(lights.empty())return nullptr;
    float currPower = 0;
    float power = u * totalPower;
    //do binary search
    for(std::size_t i = 0;i < lightPowers.size();i++){
        currPower += lightPowers[i];
        if(currPower >= power){
            return lights[i];
        }
    }
    return lights.back();
}
float PowerLightSampler::PMF(const std::shared_ptr<Light>& light) const{
    if(totalPower == 0)return 1.0f;
    return light->Power() / totalPower;
}

void PowerLightSampler::PreProcess(const AABB& bbox){
    lightPowers.clear();
    lightPowers.reserve(lights.size());
    std::vector<std::shared_ptr<Light>> validLights; //we cull lights with very small power
    validLights.reserve(lights.size());
    for(auto&& light : lights){
        light->PreProcess(bbox);
        if(light->Power() < 0.01f)continue;
        validLights.emplace_back(light);
        totalPower += lightPowers.emplace_back(light->Power());
    }
    lights = validLights;
}


