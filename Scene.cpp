#include "Scene.hpp"

bool Scene::IntersectPred(const Ray& ray, float max) const {
    return scene_bvh.IntersectPred(ray,max);
}
bool Scene::IntersectTr(Ray ray, SurfaceInteraction& interaction, glm::vec3& Tr, float max) const {
    Tr = {1,1,1};
    while(true){
        bool hit = Intersect(ray,interaction,max);
        if(ray.medium)
            Tr *= ray.medium->Tr(ray,interaction.t);
        
        if(!hit)
            return false;
        if(interaction.mat != nullptr)
            return true;
        ray = Ray(ray.at(interaction.t),ray.dir,interaction.getMedium(ray.dir));
        max-=interaction.t;
    }
    return false;
}
bool Scene::Intersect(const Ray& ray, SurfaceInteraction& interaction, float max) const {
    return scene_bvh.Intersect(ray,interaction,max);
}

void Scene::Add(const std::shared_ptr<Primitive>& ptr) {
    primitives.push_back(ptr);
}

std::vector<std::shared_ptr<Light>> Scene::GetLights() const {
    return scene_bvh.GetLights();
}

void Scene::PreProcess() {
    scene_bvh = TLAS(primitives);
}

AABB Scene::Bounding_box() const {
    return scene_bvh.Bounding_box();
}