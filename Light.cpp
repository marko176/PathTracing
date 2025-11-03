#include "Light.hpp"
#include <numbers>
#include "Sampler.hpp"
#include <thread>
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

glm::vec3 UniformInfiniteLight::L(const SurfaceInteraction& interaction, const Ray& ray) const {
    return color;
}

LightSample UniformInfiniteLight::sample(const glm::vec2& uv, float time) const{
    float z = 2.0f * uv.x - 1.0f;
    float theta = 2.0f* std::numbers::pi_v<float> * uv.y;
    float r = std::sqrt(1.0f - z*z);
    float x = r * std::cos(theta);
    float y = r * std::sin(theta);
    return {color,{{},{},SphereShape::getSphereUV(glm::normalize(glm::vec3(x,y,z)))},glm::vec3(x,y,z)};
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

glm::vec3 FunctionInfiniteLight::L(const SurfaceInteraction& interaction, const Ray& ray) const {
    return Le(ray);
}

LightSample FunctionInfiniteLight::sample(const glm::vec2& uv, float time) const{
    float z = 2.0f * uv.x - 1.0f;
    float theta = 2.0f* std::numbers::pi_v<float> * uv.y;
    float r = std::sqrt(1.0f - z*z);
    float x = r * std::cos(theta);
    float y = r * std::sin(theta);
    return {Le({{0,0,0},glm::vec3(x,y,z)}),{{},{},SphereShape::getSphereUV(glm::normalize(glm::vec3(x,y,z)))},glm::vec3(x,y,z)};
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

glm::vec3 TextureInfiniteLight::L(const SurfaceInteraction& interaction, const Ray& ray) const {
    return Le(ray);
}

LightSample TextureInfiniteLight::sample(const glm::vec2& uv, float time) const{

    float weight = random_double() * totalWeight;
    int index = std::upper_bound(accWeights.begin(), accWeights.end(), weight) - accWeights.begin();
    int cellX = index % ySamples;
    int cellY = index / ySamples;

    glm::vec2 cellUV = glm::vec2{ (cellX + uv.x) / float(xSamples),
                                  (cellY + uv.y) / float(ySamples) };
    


    float z = 2.0f * cellUV.x - 1.0f;
    float theta = 2.0f * std::numbers::pi_v<float> * cellUV.y;
    float r = std::sqrt(1.0f - z*z);
    float sx = r * std::cos(theta);
    float sy = r * std::sin(theta);
    glm::vec3 dir = glm::vec3(sx, sy, z);
    
    //float cellOmega = 4.0f * std::numbers::pi_v<float> / static_cast<float>(spp*spp);

    //float pdf = (weights[index] / totalWeight) * (1.0f / cellOmega);
    //Le({{0,0,0},dir})
    return {{0,0,0},{{},{},SphereShape::getSphereUV(glm::normalize(dir))},dir};
}

float TextureInfiniteLight::PDF(const GeometricInteraction& interaction, const Ray& ray) const{
    glm::vec3 l = Le(ray);
    float cellOmega = 4.0f * std::numbers::pi_v<float> / static_cast<float>(xSamples*ySamples);
    return (luminance(l) / totalWeight)  * (1.0f / cellOmega);
}

float TextureInfiniteLight::Power() const {
    return cachedPower;
}

void TextureInfiniteLight::PreProcess(const AABB& bbox) {
    glm::vec3 center = (bbox.max + bbox.min ) * 0.5f;
    sceneRadius = glm::distance(bbox.max,center);

    weights.clear();
    accWeights.clear();
    totalWeight = 0;
    int samples = xSamples*ySamples;
    constexpr int SPP = 64;//lower based on spp
    weights.assign(samples,0);
    accWeights.assign(samples,0);

    unsigned int threads = std::thread::hardware_concurrency();
    std::vector<std::thread> workers;
    std::atomic<int> done{0};

    auto sampler = std::make_shared<StratifiedSampler>(std::sqrt(SPP),std::sqrt(SPP));
    auto lamb = [&](){
        int k;
        std::shared_ptr<Sampler> clonedSampler = sampler->Clone();
        while((k = done.fetch_add(1,std::memory_order_relaxed))<samples){
            int y = k/ySamples;
            int x = k%ySamples;
            double temp = 0;
            for(int sp = 0;sp<SPP;sp++){
                clonedSampler->StartPixelSample({x,y},sp);
                glm::vec2 UV = clonedSampler->getPixel2D();
                glm::vec2 uv = glm::vec2{(x + UV.x) / static_cast<float>(xSamples), (y + UV.y) / static_cast<float>(ySamples)};
                float z = 2.0f * uv.x - 1.0f;
                float theta = 2.0f* std::numbers::pi_v<float> * uv.y;
                float r = std::sqrt(1.0f - z*z);
                Ray ray(glm::vec3{0,0,0},{r * std::cos(theta),r * std::sin(theta),z});
                glm::vec3 l = Le(ray);
                temp+=luminance(l);
            }
            temp/=SPP;
            weights[k]=temp;
        }
    };
    for(unsigned int t = 0;t<threads;t++){
        workers.emplace_back(lamb);
    }
    for(auto& worker : workers)worker.join();
    std::partial_sum(weights.begin(),weights.end(),accWeights.begin());
    totalWeight = accWeights.back();
    cachedPower = totalWeight / samples * powerFunction(sceneRadius);
}

bool DistantLight::isDelta() const {
    return true;
}
glm::vec3 DistantLight::L(const SurfaceInteraction& interaction, const Ray& ray) const {
    return {0,0,0};
}
LightSample DistantLight::sample(const glm::vec2& uv, float time) const{
    float z = 2.0f * uv.x - 1.0f;
    float theta = 2.0f* std::numbers::pi_v<float> * uv.y;
    float r = std::sqrt(1.0f - z*z);
    float x = r * std::cos(theta);
    float y = r * std::sin(theta);
    return {color,{{},{},uv},glm::normalize(dir + glm::vec3(x,y,z)*0.02f)};
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
glm::vec3 PointLight::L(const SurfaceInteraction& interaction, const Ray& ray) const {
    return {0,0,0};
}
LightSample PointLight::sample(const glm::vec2& uv, float time) const{
    return {color,SurfaceInteraction{p,glm::vec3{1,1,1},uv},glm::vec3{0,0,0}};
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
glm::vec3 AreaLight::L(const SurfaceInteraction& interaction, const Ray& ray) const  {
    if(oneSided && glm::dot(ray.dir,interaction.n)>0)return {0,0,0};
    return emissiveTexture->Evaluate(interaction);
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
    return cachedPower;//(color.x + color.y + color.z); FIXX!!!
}

void AreaLight::PreProcess(const AABB& bbox) {
    glm::dvec3 acc = {0,0,0};
    //have stratified sampler
    int samples = 16*16;
    StratifiedSampler sampler(std::sqrt(samples),std::sqrt(samples));
    for(int i = 0;i<samples;i++){
        acc+=emissiveTexture->Evaluate(shape->Sample(sampler.get2D()));
    }
    acc/=samples;
    cachedPower = (oneSided ? 1 : 2) * shape->Area() * luminance(acc);
}

std::shared_ptr<Shape> AreaLight::getShape() const {
    return shape;
}


//we could have transform class
//transform(interaction,ray) at some time
//light uses the transformed stuff



bool TransformedLight::isDelta() const  {
    return light->isDelta();
}
glm::vec3 TransformedLight::L(const SurfaceInteraction& interaction, const Ray& ray) const  {
    SurfaceInteraction temp;
    temp.n = normalMatrix * interaction.n;
    temp.p = transform * glm::vec4(interaction.p,1);
    return light->L(temp,ray);
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

void TransformedLight::PreProcess(const AABB& bbox) {
    light->PreProcess(bbox);
}

float TransformedLight::Power() const  {
    return light->Power();// * scale of transform;
}



bool AnimatedLight::isDelta() const  {
    return light->isDelta();
}
glm::vec3 AnimatedLight::L(const SurfaceInteraction& interaction, const Ray& ray) const  {
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

void AnimatedLight::PreProcess(const AABB& bbox) {
    light->PreProcess(bbox);
}

float AnimatedLight::Power() const  {
    return light->Power();
}