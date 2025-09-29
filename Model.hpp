#pragma once
#include "Mesh.hpp"
#include <string>
#include <vector>
#include "AABB.hpp"
#include "Primitive.hpp"
struct aiScene;
struct aiNode;
struct aiMesh;


class Material;
class Model : public Primitive{
public:
    virtual ~Model() = default;
    Model(const std::string& path);
    Model(const std::string& path,const std::shared_ptr<Material>& material, const std::shared_ptr<Medium>& medium);
    //Model(std::vector<std::shared_ptr<Mesh>> meshes);
    //buildBLAS();
    //hold full name not just path
    //fix export to have fullname.assbin
    //give virtual class for parsing to mesh
    auto GetMeshes() const -> const std::vector<std::shared_ptr<Mesh>>&;
    bool Intersect(const Ray& ray,SurfaceInteraction& interaction, float max) const override{
        return model_bvh.Intersect(ray,interaction,max);
    }
    AABB BoundingBox() const override{
        return model_bvh.BoundingBox();
    }
    
    bool IntersectPred(const Ray& ray, float max) const override{
        return model_bvh.IntersectPred(ray,max);
    }
    
    std::vector<std::shared_ptr<Light>> GetLights() const override {
        return model_bvh.GetLights();
    }

    
private:
    auto load_model(const std::string& path) -> bool;
    auto process_node(aiNode* node, const aiScene* scene) -> void;
    auto process_mesh(aiMesh* mesh, const aiScene* scene) -> std::shared_ptr<Mesh>;

    std::vector<std::shared_ptr<Mesh>> meshes;
    BLAS model_bvh;
    std::string model_path;
};