#include "Primitive.hpp"

AABB GeometricPrimitive::BoundingBox() const{
    return shape->BoundingBox();
}
bool GeometricPrimitive::IntersectPred(const Ray& ray, float max) const{
    if(material && material->HasAlpha()){
        SurfaceInteraction intr;
        bool hit = shape->Intersect(ray, intr, max);
        return hit && material->Alpha(intr.uv);
    } else{
        return shape->IntersectPred(ray, max);
    }
}
bool GeometricPrimitive::Intersect(const Ray& ray, SurfaceInteraction& interaction, float max) const{
    SurfaceInteraction intr;
    bool hit = shape->Intersect(ray, intr, max);
    if(!hit || (material && !material->Alpha(intr.uv)))return false;
    
    interaction = intr;
    interaction.AreaLight = areaLight;
    interaction.mat = material;
    interaction.medium = medium;

    return true;
}

std::vector<std::shared_ptr<Light>> GeometricPrimitive::GetLights() const{
    return areaLight != nullptr ? std::vector<std::shared_ptr<Light>>{areaLight} : std::vector<std::shared_ptr<Light>> {};
}

AABB TransformedPrimitive::BoundingBox() const{
    AABB bbox;
    AABB temp = primitive->BoundingBox();

    for(int i = 0;i < 8;i++){
        bbox.Expand(transform * glm::vec4(temp.Corner(i), 1));
    }
    return bbox;
}
bool TransformedPrimitive::IntersectPred(const Ray& ray, float max) const{
    glm::vec3 dir = invTransform * glm::vec4(ray.dir, 0);
    float length = glm::length(dir);
    Ray transformedRay(invTransform * glm::vec4(ray.origin, 1), dir / length, ray.time);
    return primitive->IntersectPred(transformedRay, max * length);
}
bool TransformedPrimitive::Intersect(const Ray& ray, SurfaceInteraction& interaction, float max) const{
    glm::vec3 dir = invTransform * glm::vec4(ray.dir, 0);
    float length = glm::length(dir);
    Ray transformedRay(invTransform * glm::vec4(ray.origin, 1), dir / length, ray.time);
    SurfaceInteraction temp;
    if(!primitive->Intersect(transformedRay, temp, max * length))return false;

    glm::mat3 normalMatrix = glm::transpose(glm::inverse(glm::mat3(transform)));
    interaction.p = transform * glm::vec4(temp.p, 1);
    interaction.AreaLight = temp.AreaLight;
    interaction.mat = temp.mat;
    interaction.n = glm::normalize(normalMatrix * temp.n);
    interaction.ns = glm::normalize(normalMatrix * temp.ns);
    interaction.t = temp.t / length;
    interaction.tangent = glm::normalize(transform * glm::vec4(temp.tangent, 0));
    interaction.uv = temp.uv;
    interaction.medium = temp.medium;

    return true;
}
std::vector<std::shared_ptr<Light>> TransformedPrimitive::GetLights() const{
    std::vector<std::shared_ptr<Light>> temp;
    for(const std::shared_ptr<Light>& l : primitive->GetLights()){
        temp.push_back(std::make_shared<TransformedLight>(l, transform));
    }
    return temp;
}



AABB AnimatedPrimitive::BoundingBox() const{
    AABB bbox = primitive->BoundingBox();
    bbox.Expand(TransformedPrimitive(primitive, glm::translate(glm::mat4(1), dir)).BoundingBox());
    return bbox;
}
bool AnimatedPrimitive::IntersectPred(const Ray& ray, float max) const{
    float t = glm::clamp(ray.time - timeBounds.x, timeBounds.x, timeBounds.y) / (timeBounds.y - timeBounds.x);
    return TransformedPrimitive(primitive, glm::translate(glm::mat4(1), dir * t)).IntersectPred(ray, max);
}
bool AnimatedPrimitive::Intersect(const Ray& ray, SurfaceInteraction& interaction, float max) const{
    float t = glm::clamp(ray.time - timeBounds.x, timeBounds.x, timeBounds.y) / (timeBounds.y - timeBounds.x);
    return TransformedPrimitive(primitive, glm::translate(glm::mat4(1), dir * t)).Intersect(ray, interaction, max);
}
std::vector<std::shared_ptr<Light>> AnimatedPrimitive::GetLights() const{
    std::vector<std::shared_ptr<Light>> temp;
    for(const std::shared_ptr<Light>& l : primitive->GetLights()){
        temp.push_back(std::make_shared<AnimatedLight>(l, dir, timeBounds));
    }
    return temp;
}