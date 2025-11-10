#pragma once
#include "AABB.hpp"
#include "Primitive.hpp"
#include "BVH.hpp"
class Medium;
class Scene{
public:
    Scene() = default;
    Scene(const std::shared_ptr<Medium>& medium) : sceneMedium(medium){}

    bool IntersectPred(const Ray& ray, float max = std::numeric_limits<float>::infinity()) const;

    bool IntersectTr(Ray ray, SurfaceInteraction& interaction, glm::vec3& Tr, float max = std::numeric_limits<float>::infinity()) const;

    bool Intersect(const Ray& ray, SurfaceInteraction& interaction, float max = std::numeric_limits<float>::infinity()) const;

    void Add(const std::shared_ptr<Primitive>& ptr);

    std::vector<std::shared_ptr<Light>> GetLights() const;

    template <typename T>
    void BuildTlas(){
        scene_bvh = std::make_shared<T>(std::move(primitives));
    }

    std::shared_ptr<Medium> GetMedium() const;

    void SetMedium(const std::shared_ptr<Medium>& medium);

    AABB BoundingBox() const;

    std::vector<std::shared_ptr<InfiniteLight>> infiniteLights;
private:
    std::shared_ptr<BVHBase<std::shared_ptr<Primitive>>> scene_bvh;//std::shared_ptr<BvhBase<std::shared_ptr<Primitive>> 
    std::vector<std::shared_ptr<Primitive>> primitives;
    std::shared_ptr<Medium> sceneMedium;
};