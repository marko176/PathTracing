#pragma once
#include "Hit_record.hpp"
#include "Primitive.hpp"
class Scene {
public:

    bool IntersectPred(const Ray& ray, float max = std::numeric_limits<float>::infinity()) const;

    bool IntersectTr(Ray ray, SurfaceInteraction& interaction, glm::vec3& Tr, float max = std::numeric_limits<float>::infinity()) const;

    bool Intersect(const Ray& ray, SurfaceInteraction& interaction, float max = std::numeric_limits<float>::infinity()) const;

    void Add(const std::shared_ptr<Primitive>& ptr);
    
    std::vector<std::shared_ptr<Light>> GetLights() const;
    
    void PreProcess();
    
    AABB BoundingBox() const;

    TLAS scene_bvh;//tlas bvh should take shared ptrs / unique ptrs and own them
private:
    std::vector<std::shared_ptr<Primitive>> primitives;
};