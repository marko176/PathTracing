#include "Scene.hpp"
#include "Medium.hpp"
bool Scene::IntersectPred(const Ray& ray, float max) const {
    return scene_bvh.IntersectPred(ray,max);
}
bool Scene::IntersectTr(Ray ray, SurfaceInteraction& interaction, glm::vec3& Tr, float max) const {
    Tr = {1,1,1};
    while(max > 0){
        bool hit = Intersect(ray,interaction,max);

        if(!hit){
            if(ray.medium)
                Tr *= ray.medium->Tr(ray,max);
                
            return false;
        }

        if(ray.medium)
            Tr *= ray.medium->Tr(ray,interaction.t);

        if(interaction.mat != nullptr)
            return true;
        ray = Ray(ray.at(interaction.t),ray.dir,ray.time,interaction.getMedium(ray.dir));
        max-= interaction.t;
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

    std::vector<std::shared_ptr<Light>> lights = scene_bvh.GetLights();
    lights.insert(lights.end(),infiniteLights.begin(),infiniteLights.end());
    return lights; //if we want to sample the skybox
    return scene_bvh.GetLights();
}

void Scene::PreProcess() {
    scene_bvh = TLAS(primitives);
}

AABB Scene::BoundingBox() const {
    return scene_bvh.BoundingBox();
}

std::shared_ptr<Medium> Scene::GetMedium() const {
    return sceneMedium;
}

void Scene::SetMedium(const std::shared_ptr<Medium>& medium) {
    sceneMedium = medium;
}