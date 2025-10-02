#include "Light.hpp"
#include <numbers>
#include "Sampler.hpp"

inline double luminance(const glm::dvec3& v){//taken from gilm to film
    return dot(v, glm::dvec3(0.2126, 0.7152, 0.0722));
}


void InfiniteLight::PreProcess(const AABB& bbox) {
    glm::vec3 center = (bbox.max + bbox.min ) * 0.5f;
    sceneRadius = glm::distance(bbox.max,center);
}

bool InfiniteLight::isDelta() const {
    return false;
}

float InfiniteLight::PDF(const GeometricInteraction& interaction, float time) const{
    //1/area => 1/inf = 0
    return 0;
}

float InfiniteLight::PDF(const GeometricInteraction& interaction, const Ray& ray) const{
    return 0;
}

glm::vec3 UniformInfiniteLight::Le(const Ray& ray) const {
    return color;
}

glm::vec3 UniformInfiniteLight::L(const GeometricInteraction& interaction, const Ray& ray) const {
    return color;
}

LightSample UniformInfiniteLight::sample(const glm::vec2& uv, float time) const{
    float z = 2.0f * uv.x - 1.0f;
    float theta = 2.0f* std::numbers::pi_v<float> * uv.y;
    float r = std::sqrt(1.0f - z*z);
    float x = r * std::cos(theta);
    float y = r * std::sin(theta);
    return {color,{},glm::vec3(x,y,z)};
}

float UniformInfiniteLight::PDF(const GeometricInteraction& interaction, const Ray& ray) const{
    return 1.0f / (4.0f * std::numbers::pi_v<float>);
}

float UniformInfiniteLight::Power() const {
    return (color.x + color.y + color.z) * powerFunction(sceneRadius);//
}


glm::vec3 FunctionInfiniteLight::Le(const Ray& ray) const {
    return lightFunction(ray);
}

glm::vec3 FunctionInfiniteLight::L(const GeometricInteraction& interaction, const Ray& ray) const {
    return Le(ray);
}

LightSample FunctionInfiniteLight::sample(const glm::vec2& uv, float time) const{
    float z = 2.0f * uv.x - 1.0f;
    float theta = 2.0f* std::numbers::pi_v<float> * uv.y;
    float r = std::sqrt(1.0f - z*z);
    float x = r * std::cos(theta);
    float y = r * std::sin(theta);
    return {Le({{0,0,0},glm::vec3(x,y,z)}),{},glm::vec3(x,y,z)};
}

float FunctionInfiniteLight::PDF(const GeometricInteraction& interaction, const Ray& ray) const{
    //wrong if we sample based on funciton
    return 1.0f / (4.0f * std::numbers::pi_v<float>);
}

float FunctionInfiniteLight::Power() const {
    return cachedPower;
}

void FunctionInfiniteLight::PreProcess(const AABB& bbox) {
    glm::vec3 center = (bbox.max + bbox.min ) * 0.5f;
    sceneRadius = glm::distance(bbox.max,center);

    double acc = 0;
    int samples = 100*100;
    int sqrtSamples = std::sqrt(samples);
    constexpr int SPP = 100;//lower based on spp
    auto sampler = std::make_shared<StratifiedSampler>(sqrtSamples,sqrtSamples);
    for(int y = 0;y < sqrtSamples; y++){
        for(int x = 0; x < sqrtSamples; x++){
            double temp = 0;
            sampler->StartPixelSample({x,y},x + y*sqrtSamples);
            for(int sp = 0;sp<SPP;sp++){
                glm::vec2 UV = sampler->get2D();
                glm::vec2 uv = glm::vec2{(x + UV.x), (y + UV.y)} / static_cast<float>(sqrtSamples);
                float z = 2.0f * uv.x - 1.0f;
                float theta = 2.0f* std::numbers::pi_v<float> * uv.y;
                float r = std::sqrt(1.0f - z*z);
                Ray ray(glm::vec3{0,0,0},{r * std::cos(theta),r * std::sin(theta),z});
                glm::vec3 l = Le(ray);
                temp+=luminance(l);
            }
            temp/=SPP;
            acc += temp;
        }
    }
    cachedPower = acc / samples * powerFunction(sceneRadius);
}


glm::vec3 TextureInfiniteLight::Le(const Ray& ray) const {
    return LeScale*tex->Evaluate({{0,0,0},{0,0,0},SphereShape::getSphereUV(ray.dir)});//must be normalized
}

glm::vec3 TextureInfiniteLight::L(const GeometricInteraction& interaction, const Ray& ray) const {
    return Le(ray);
}

LightSample TextureInfiniteLight::sample(const glm::vec2& uv, float time) const{

    float weight = random_double() * totalWeight;
    int index = std::upper_bound(accWeights.begin(), accWeights.end(), weight) - accWeights.begin();
    int cellX = index % spp;
    int cellY = index / spp;

    glm::vec2 cellUV = glm::vec2{ (cellX + uv.x) / float(spp),
                                  (cellY + uv.y) / float(spp) };
    


    float z = 2.0f * cellUV.x - 1.0f;
    float theta = 2.0f * std::numbers::pi_v<float> * cellUV.y;
    float r = std::sqrt(1.0f - z*z);
    float sx = r * std::cos(theta);
    float sy = r * std::sin(theta);
    glm::vec3 dir = glm::vec3(sx, sy, z);
    
    //float cellOmega = 4.0f * std::numbers::pi_v<float> / static_cast<float>(spp*spp);

    //float pdf = (weights[index] / totalWeight) * (1.0f / cellOmega);
    //Le({{0,0,0},dir})
    return {{0,0,0},{},dir};
}

float TextureInfiniteLight::PDF(const GeometricInteraction& interaction, const Ray& ray) const{
    glm::vec3 l = Le(ray);
    float cellOmega = 4.0f * std::numbers::pi_v<float> / static_cast<float>(spp*spp);
    return (luminance(l) / totalWeight)  * (1.0f / cellOmega);
}

float TextureInfiniteLight::Power() const {
    return cachedPower;
}

void TextureInfiniteLight::PreProcess(const AABB& bbox) {
    glm::vec3 center = (bbox.max + bbox.min ) * 0.5f;
    sceneRadius = glm::distance(bbox.max,center);

    weights.clear();
    totalWeight = 0;
    double acc = 0;
    int samples = spp*spp;
    int sqrtSamples = std::sqrt(samples);
    constexpr int SPP = 100;//lower based on spp
    auto sampler = std::make_shared<StratifiedSampler>(sqrtSamples,sqrtSamples);
    for(int y = 0;y < sqrtSamples; y++){
        for(int x = 0; x < sqrtSamples; x++){
            double temp = 0;
            sampler->StartPixelSample({x,y},x + y*sqrtSamples);
            for(int sp = 0;sp<SPP;sp++){
                glm::vec2 UV = sampler->get2D();
                glm::vec2 uv = glm::vec2{(x + UV.x), (y + UV.y)} / static_cast<float>(sqrtSamples);
                float z = 2.0f * uv.x - 1.0f;
                float theta = 2.0f* std::numbers::pi_v<float> * uv.y;
                float r = std::sqrt(1.0f - z*z);
                Ray ray(glm::vec3{0,0,0},{r * std::cos(theta),r * std::sin(theta),z});
                glm::vec3 l = Le(ray);
                temp+=luminance(l);
            }
            temp/=SPP;
            acc += temp;
            weights.push_back(temp);
            totalWeight += weights.back();
            accWeights.push_back(totalWeight);
        }
    }
    cachedPower = acc / samples * powerFunction(sceneRadius);
    
}

bool DistantLight::isDelta() const {
    return true;
}
glm::vec3 DistantLight::L(const GeometricInteraction& interaction, const Ray& ray) const {
    return {0,0,0};
}
LightSample DistantLight::sample(const glm::vec2& uv, float time) const{
    float z = 2.0f * uv.x - 1.0f;
    float theta = 2.0f* std::numbers::pi_v<float> * uv.y;
    float r = std::sqrt(1.0f - z*z);
    float x = r * std::cos(theta);
    float y = r * std::sin(theta);
    return {color,{},glm::normalize(dir + glm::vec3(x,y,z)*0.02f)};
}
float DistantLight::PDF(const GeometricInteraction& interaction, float time) const{
    return 0;
}
float DistantLight::PDF(const GeometricInteraction& interaction, const Ray& ray) const{
    return 0;
}
float DistantLight::Power() const {
    return (color.x + color.y + color.z) * powerFunction(sceneRadius);//
}
void DistantLight::PreProcess(const AABB& bbox) {
    glm::vec3 center = (bbox.max + bbox.min ) * 0.5f;
    sceneRadius = glm::distance(bbox.max,center);
}

bool PointLight::isDelta() const {
    return true;
}
glm::vec3 PointLight::L(const GeometricInteraction& interaction, const Ray& ray) const {
    return {0,0,0};
}
LightSample PointLight::sample(const glm::vec2& uv, float time) const{
    return {color,GeometricInteraction{p,glm::vec3{1,1,1}},glm::vec3{0,0,0}};
}
float PointLight::PDF(const GeometricInteraction& interaction, float time) const{
    return 0;
}
float PointLight::PDF(const GeometricInteraction& interaction, const Ray& ray) const{
    return 0;
}
float PointLight::Power() const {
    return (color.x + color.y + color.z) * powerFunction(sceneRadius);//
}
void PointLight::PreProcess(const AABB& bbox) {
    glm::vec3 center = (bbox.max + bbox.min ) * 0.5f;
    sceneRadius = glm::distance(bbox.max,center);
}


bool AreaLight::isDelta() const  {
    return false;
}
glm::vec3 AreaLight::L(const GeometricInteraction& interaction, const Ray& ray) const  {
    if(oneSided && glm::dot(ray.dir,interaction.n)>0)return {0,0,0};
    return color;
}
LightSample AreaLight::sample(const glm::vec2& uv, float time) const {
    return {{0,0,0},shape->Sample(uv),{}};
}
float AreaLight::PDF(const GeometricInteraction& interaction, float time) const {
    return shape->PDF(interaction);//this is wrong we should co dist^2/cosine...
}
float AreaLight::PDF(const GeometricInteraction& interaction, const Ray& ray) const {
    if(oneSided){
        return glm::dot(-ray.dir,interaction.n) > 0 ? shape->PDF(interaction,ray) : 0;
    }
    return shape->PDF(interaction,ray);
}
float AreaLight::Power() const  {
    return (oneSided ? 1 : 2) * shape->Area() * (color.x + color.y + color.z);
}

std::shared_ptr<Shape> AreaLight::getShape() const {
    return shape;
}


//we could ahve transform class
//transform(interaction,ray) at some time
//ligth uses the transformed stuff



bool TransformedLight::isDelta() const  {
    return light->isDelta();
}
glm::vec3 TransformedLight::L(const GeometricInteraction& interaction, const Ray& ray) const  {
    GeometricInteraction temp;
    temp.n = normalMatrix * interaction.n;
    temp.p = transform * glm::vec4(interaction.p,1);
    return light->L(interaction,ray);
}
LightSample TransformedLight::sample(const glm::vec2& uv, float time) const {
    LightSample lightSample = light->sample(uv,time);
    lightSample.dir = transform * glm::vec4(lightSample.dir,0);
    lightSample.interaction.p = transform * glm::vec4(lightSample.interaction.p,1);
    lightSample.interaction.n = normalMatrix * lightSample.interaction.n;
    return lightSample;
}
float TransformedLight::PDF(const GeometricInteraction& interaction, float time) const {
    glm::vec3 p = glm::vec3(invTransform * glm::vec4(interaction.p, 1.0f));
    glm::vec3 n = glm::normalize(glm::vec3(invTransform * glm::vec4(interaction.n, 0.0f)));

    return light->PDF(GeometricInteraction(p,n),time);
}
float TransformedLight::PDF(const GeometricInteraction& interaction, const Ray& ray) const {
    glm::vec3 p = glm::vec3(invTransform * glm::vec4(interaction.p, 1.0f));
    glm::vec3 n = glm::normalize(glm::vec3(invTransform * glm::vec4(interaction.n, 0.0f)));
    Ray localRay(glm::vec3(invTransform * glm::vec4(ray.origin,1)),
                    glm::normalize(glm::vec3(invTransform * glm::vec4(ray.dir,0))),ray.time);

    return light->PDF(GeometricInteraction(p,n),localRay);
}
float TransformedLight::Power() const  {
    return light->Power();
}



bool AnimatedLight::isDelta() const  {
    return light->isDelta();
}
glm::vec3 AnimatedLight::L(const GeometricInteraction& interaction, const Ray& ray) const  {
    float t = glm::clamp(ray.time - timeBounds.x,timeBounds.x,timeBounds.y) / (timeBounds.y - timeBounds.x);
    return TransformedLight(light,glm::translate(glm::mat4(1),dir * t)).L(interaction,ray);
}
LightSample AnimatedLight::sample(const glm::vec2& uv, float time) const {
    float t = glm::clamp(time - timeBounds.x,timeBounds.x,timeBounds.y) / (timeBounds.y - timeBounds.x);
    return TransformedLight(light,glm::translate(glm::mat4(1),dir * t)).sample(uv,time);
}
float AnimatedLight::PDF(const GeometricInteraction& interaction, float time) const {
    float t = glm::clamp(time - timeBounds.x,timeBounds.x,timeBounds.y) / (timeBounds.y - timeBounds.x);
    return TransformedLight(light,glm::translate(glm::mat4(1),dir * t)).PDF(interaction,time);
}
float AnimatedLight::PDF(const GeometricInteraction& interaction, const Ray& ray) const {
    float t = glm::clamp(ray.time - timeBounds.x,timeBounds.x,timeBounds.y) / (timeBounds.y - timeBounds.x);
    return TransformedLight(light,glm::translate(glm::mat4(1),dir * t)).PDF(interaction,ray);
}
float AnimatedLight::Power() const  {
    return light->Power();
}