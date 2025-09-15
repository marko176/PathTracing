#pragma once
#include "Mesh.hpp"
#include <string>
#include <vector>
#include "Hit_record.hpp"
#include "Primitive.hpp"
struct aiScene;
struct aiNode;
struct aiMesh;



class Model : public Primitive{
    public:
    Model(const std::string& path);
    auto get_meshes() const -> const std::vector<Mesh>&;
    bool Intersect(const Ray& ray,SurfaceInteraction& interaction, float max) const override{
        return model_bvh.Intersect(ray,interaction,max);
    }
    AABB Bounding_box() const override{
        return model_bvh.Bounding_box();
    }
    
    bool IntersectPred(const Ray& ray, float max) const override{
        return model_bvh.IntersectPred(ray,max);
    }
    
    std::vector<std::shared_ptr<Light>> GetLights() const override {
        return model_bvh.GetLights();
    }
    private:
    auto load_model(std::string path) -> bool;
    auto process_node(aiNode* node, const aiScene* scene) -> void;
    auto process_mesh(aiMesh* mesh, const aiScene* scene) -> Mesh;
    //std::vector<BVH> mesh_bvhs;
    //std::vector<GeometricPrimitive> primitives;
    std::vector<Mesh> meshes;// move to private!!!!
    BLAS model_bvh;//was primitiveBVH
    std::string model_path;


    //std::shared_ptr<int> tempPTR;
    //TLAS_BVH_Prim model_bvh;

};