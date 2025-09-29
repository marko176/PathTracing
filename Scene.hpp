#pragma once
#include "AABB.hpp"
#include "Primitive.hpp"
class Medium;
class Scene {
public:
    Scene() = default;
    Scene(const std::shared_ptr<Medium>& medium) : sceneMedium(medium) {}

    bool IntersectPred(const Ray& ray, float max = std::numeric_limits<float>::infinity()) const;

    bool IntersectTr(Ray ray, SurfaceInteraction& interaction, glm::vec3& Tr, float max = std::numeric_limits<float>::infinity()) const;

    bool Intersect(const Ray& ray, SurfaceInteraction& interaction, float max = std::numeric_limits<float>::infinity()) const;

    void Add(const std::shared_ptr<Primitive>& ptr);
    
    std::vector<std::shared_ptr<Light>> GetLights() const;
    
    void PreProcess();

    std::shared_ptr<Medium> GetMedium() const;

    void SetMedium(const std::shared_ptr<Medium>& medium);

    AABB BoundingBox() const;

private:
    TLAS scene_bvh;//tlas bvh should take shared ptrs / unique ptrs and own them
    std::vector<std::shared_ptr<Primitive>> primitives;
    std::shared_ptr<Medium> sceneMedium;
};