#include "Primitive.hpp"
#include <chrono>
AABB GeometricPrimitive::BoundingBox() const {
    return shape->BoundingBox();
}
bool GeometricPrimitive::IntersectPred(const Ray& ray, float max) const {
    bool hit = shape->IntersectPred(ray,max);
    //if(material->Alpha()) we need UV coords from intersecpt pred!!
    return hit;
}
bool GeometricPrimitive::Intersect(const Ray& ray, SurfaceInteraction& interaction, float max) const {
    bool hit = shape->Intersect(ray,interaction,max);
    if(hit){
        interaction.AreaLight = areaLight;
        interaction.mat = material;
        interaction.medium = medium;
    }
    return hit;
}

std::vector<std::shared_ptr<Light>> GeometricPrimitive::GetLights() const {
    return areaLight != nullptr ? std::vector<std::shared_ptr<Light>>{areaLight} : std::vector<std::shared_ptr<Light>>{};
}

AABB TransformedPrimitive::BoundingBox() const {
    AABB bbox;
    AABB temp = primitive->BoundingBox();
    for(int i = 0;i<8;i++){
        bbox.Expand(transform * glm::vec4(temp.Corner(i),1));
    }
    return bbox;
}
bool TransformedPrimitive::IntersectPred(const Ray& ray, float max) const {
    glm::vec3 dir = invTransform * glm::vec4(ray.dir,0);
    float length = glm::length(dir);
    Ray transformedRay(invTransform * glm::vec4(ray.origin,1), dir / length);
    return primitive->IntersectPred(transformedRay,max * length);
}
bool TransformedPrimitive::Intersect(const Ray& ray, SurfaceInteraction& interaction, float max) const {
    glm::vec3 dir = invTransform * glm::vec4(ray.dir,0);
    float length = glm::length(dir);
    Ray transformedRay(invTransform * glm::vec4(ray.origin,1), dir / length);
    SurfaceInteraction temp;
    if(!primitive->Intersect(transformedRay,temp,max * length))return false;

    glm::mat3 normalMatrix = glm::transpose(glm::inverse(glm::mat3(transform)));
    interaction.p = transform * glm::vec4(temp.p,1);
    interaction.AreaLight = temp.AreaLight;
    interaction.bitangent = transform * glm::vec4(temp.bitangent,0);
    if(interaction.bitangent != glm::vec3(0,0,0))interaction.bitangent = glm::normalize(interaction.bitangent);
    interaction.mat = temp.mat;
    interaction.n = glm::normalize(normalMatrix * temp.n);
    interaction.ns = glm::normalize(normalMatrix * temp.ns);
    interaction.t = temp.t / length;
    interaction.tangent = transform * glm::vec4(temp.tangent,0);
    if(interaction.tangent != glm::vec3(0,0,0))interaction.tangent = glm::normalize(interaction.tangent);
    interaction.uv = temp.uv;
    interaction.medium = temp.medium;

    return true;
}
std::vector<std::shared_ptr<Light>> TransformedPrimitive::GetLights() const {
    std::vector<std::shared_ptr<Light>> temp;
    for(const std::shared_ptr<Light>& l : primitive->GetLights()){
        temp.push_back(std::make_shared<TransformedLight>(transform,l));
    }
    return temp;
}
