#include "Light.hpp"

bool InfiniteLight::isDelta() const {
    return true;
}
glm::vec3 InfiniteLight::L(const GeometricInteraction& interaction, const Ray& ray) const {
    return color;
}
LightSample InfiniteLight::sample(const glm::vec2& uv) const{
    float z = 2.0f * uv.x - 1.0f;
    float theta = 2.0f* std::numbers::pi_v<float> * uv.y;
    float r = std::sqrt(1.0f - z*z);
    float x = r * std::cos(theta);
    float y = r * std::sin(theta);
    return {{},glm::normalize(dir + glm::vec3(x,y,z)*0.02f)};
}
float InfiniteLight::PDF(const GeometricInteraction& interaction) const{
    return 0;
}
float InfiniteLight::PDF(const GeometricInteraction& interaction, const Ray& ray) const{
    return 0;
}
float InfiniteLight::Power() const {
    return (color.x + color.y + color.z) * powerFunction(sceneRadius);//
}
void InfiniteLight::PreProcess(const AABB& bbox) {
    glm::vec3 center = (bbox.max + bbox.min ) * 0.5f;
    sceneRadius = glm::distance(bbox.max,center);
}

bool PointLight::isDelta() const {
    return true;
}
glm::vec3 PointLight::L(const GeometricInteraction& interaction, const Ray& ray) const {
    return color;
}
LightSample PointLight::sample(const glm::vec2& uv) const{
    return {GeometricInteraction{p,glm::vec3{1,1,1}},glm::vec3{0,0,0}};
}
float PointLight::PDF(const GeometricInteraction& interaction) const{
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
LightSample AreaLight::sample(const glm::vec2& uv) const {
    return {shape->Sample(uv),{}};
}
float AreaLight::PDF(const GeometricInteraction& interaction) const {
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


bool TransformedLight::isDelta() const  {
    return light->isDelta();
}
glm::vec3 TransformedLight::L(const GeometricInteraction& interaction, const Ray& ray) const  {
    GeometricInteraction temp;
    temp.n = normalMatrix * interaction.n;
    temp.p = transform * glm::vec4(interaction.p,1);
    return light->L(interaction,ray);
}
LightSample TransformedLight::sample(const glm::vec2& uv) const {
    LightSample lightSample = light->sample(uv);
    lightSample.dir = transform * glm::vec4(lightSample.dir,0);
    lightSample.interaction.p = transform * glm::vec4(lightSample.interaction.p,1);
    lightSample.interaction.n = normalMatrix * lightSample.interaction.n;
    return lightSample;
}
float TransformedLight::PDF(const GeometricInteraction& interaction) const {
    glm::vec3 p = glm::vec3(invTransform * glm::vec4(interaction.p, 1.0f));
    glm::vec3 n = glm::normalize(glm::vec3(invTransform * glm::vec4(interaction.n, 0.0f)));

    // Compute object-space PDF
    return light->PDF(GeometricInteraction(p,n));
}
float TransformedLight::PDF(const GeometricInteraction& interaction, const Ray& ray) const {
    glm::vec3 p = glm::vec3(invTransform * glm::vec4(interaction.p, 1.0f));
    glm::vec3 n = glm::normalize(glm::vec3(invTransform * glm::vec4(interaction.n, 0.0f)));
    Ray localRay(glm::vec3(invTransform * glm::vec4(ray.origin,1)),
                    glm::normalize(glm::vec3(invTransform * glm::vec4(ray.dir,0))));

    return light->PDF(GeometricInteraction(p,n),localRay);
}
float TransformedLight::Power() const  {
    return light->Power();
}