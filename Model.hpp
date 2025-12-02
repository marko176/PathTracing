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
    //Model(std::vector<std::shared_ptr<Mesh>> meshes);
    //buildBLAS();
    //hold full name not just path
    //fix export to have fullname.assbin
    //give virtual class for parsing to mesh
    std::vector<std::shared_ptr<Mesh>> GetMeshes() const;

    bool Intersect(const Ray& ray, SurfaceInteraction& interaction, float max) const override{
        return model_bvh->Intersect(ray, interaction, max);
    }

    bool IntersectPred(const Ray& ray, float max) const override{
        return model_bvh->IntersectPred(ray, max);
    }

    AABB BoundingBox() const override{
        return model_bvh->BoundingBox();
    }

    std::vector<std::shared_ptr<Light>> GetLights() const override{
        return model_bvh->GetLights();
    }

    template <typename T>
        requires std::is_base_of_v<BVHBase<GeometricPrimitive>, T>
    void BuildBlas(){
        std::vector<GeometricPrimitive> primitives;
        primitives.reserve(1'000'000);

        for(const std::shared_ptr<Mesh>& m : meshes){
            uint32_t n = m->triangle_count;
            for(uint32_t j = 0;j < n;j++){
                std::shared_ptr<Shape> shape = std::shared_ptr<Shape>(m->GetControlPtr(), m->GetShape(j));
                std::shared_ptr<AreaLight> area = m->GetEmissiveTexture() != nullptr ? std::make_shared<AreaLight>(shape, m->GetEmissiveTexture()) : nullptr;
                if(area){
                    area->PreProcess({});
                    if(area->Power() <= std::numeric_limits<float>::epsilon())area = nullptr;
                }
                primitives.emplace_back(shape, m->GetMaterial(), area, m->GetMedium());
            }
        }
        model_bvh = std::make_shared<T>(std::move(primitives));
    }

    template <typename T>
        requires std::is_base_of_v<BVHBase<GeometricPrimitive>, T>
    void BuildBlas(const std::shared_ptr<Material>& material, const std::shared_ptr<Medium>& medium){
        std::vector<GeometricPrimitive> primitives;
        primitives.reserve(1'000'000);

        for(const std::shared_ptr<Mesh>& m : meshes){
            uint32_t n = m->triangle_count;
            for(uint32_t j = 0;j < n;j++){
                std::shared_ptr<Shape> shape = std::shared_ptr<Shape>(m->GetControlPtr(), m->GetShape(j));
                std::shared_ptr<AreaLight> area = m->GetEmissiveTexture() != nullptr ? std::make_shared<AreaLight>(shape, m->GetEmissiveTexture()) : nullptr;
                if(area){
                    area->PreProcess({});
                    if(area->Power() <= std::numeric_limits<float>::epsilon())area = nullptr;
                }
                primitives.emplace_back(shape, material, area, medium);
            }
        }
        model_bvh = std::make_shared<T>(std::move(primitives));
    }
private:


    std::vector<std::shared_ptr<Texture>> GetTextures(aiMaterial* material) const;

    //std::shared_ptr<Material> SetupMaterial(aiMaterial* material) const;

    std::shared_ptr<Material> SetupMaterial(const std::vector<std::shared_ptr<Texture>>& textures, aiMaterial* material) const;

    static std::string GetFormat(const std::string& filepath){
        std::size_t index = filepath.find_last_of('.');
        if(index == std::string::npos)return "";
        std::string tmp = filepath.substr(index + 1);
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), [](unsigned char c){
            return std::tolower(c);
            });
        return tmp;
    }

    static std::string GetModelDirectory(const std::string& filepath){
        return filepath.substr(0, filepath.find_last_of('/')).append("/");
    }

    static std::string GetModelPath(const std::string& filepath){
        std::size_t index = filepath.find_last_of('.');
        if(index == std::string::npos)return filepath;
        return filepath.substr(0, index);
    }

    bool load_model(const std::string& path);
    void process_node(aiNode* node, const aiScene* scene);
    std::shared_ptr<Mesh> process_mesh(aiMesh* mesh, const aiScene* scene);

    std::vector<std::shared_ptr<Mesh>> meshes;
    std::shared_ptr<BVHBase<GeometricPrimitive>> model_bvh;//shared_ptr<bvh_base<T>>
    std::string model_path;
};