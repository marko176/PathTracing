#pragma once
#include "Mesh.hpp"
#include <string>
#include <vector>
#include "AABB.hpp"
#include "Primitive.hpp"
#include "BVH.hpp"
struct aiScene;
struct aiNode;
struct aiMesh;
class aiMaterial;





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


    std::vector<std::shared_ptr<Texture>> GetTextures(aiMaterial* material) const;

    //std::shared_ptr<Material> SetupMaterial(aiMaterial* material) const;

    std::shared_ptr<Material> SetupOBJMaterial(const std::vector<std::shared_ptr<Texture>>& textures,aiMaterial* material) const;

    static std::string GetFormat(const std::string& filepath){
        std::size_t index = filepath.find_last_of('.');
        if(index == std::string::npos)return "";
        std::string tmp = filepath.substr(index+1);
        std::transform(tmp.begin(),tmp.end(),tmp.begin(),[](unsigned char c){
            return std::tolower(c);
        });
        return tmp;
    }

    static std::string GetModelDirectory(const std::string& filepath){
        return filepath.substr(0,filepath.find_last_of('/')).append("/");
    }

    static std::string GetModelPath(const std::string& filepath){
        std::size_t index = filepath.find_last_of('.');
        if(index == std::string::npos)return filepath;
        return filepath.substr(0,index);
    }

    bool load_model(const std::string& path);
    void process_node(aiNode* node, const aiScene* scene);
    std::shared_ptr<Mesh> process_mesh(aiMesh* mesh, const aiScene* scene);

    std::vector<std::shared_ptr<Mesh>> meshes;
    BLAS4 model_bvh;//shared_ptr<bvh_base<T>>
    std::string model_path;
};