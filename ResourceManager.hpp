#pragma once
#include <string>
#include <unordered_map>
#include <cstdint>
#include "Texture.hpp"
#include "Model.hpp"
#include <memory>
#include <filesystem>
class ResourceManager{
public:
    static ResourceManager& get_instance(){
        static ResourceManager instance;
        return instance;
    }

    ResourceManager(const ResourceManager&) = delete;
    ResourceManager& operator=(const ResourceManager&) = delete;

    [[nodiscard]] std::shared_ptr<Texture> GetImageTexture(const std::string& path, bool gammaCorrection = false, const glm::vec3& colorScale = glm::vec3(1), bool invert = false);//maybe get image?


    template <typename T, typename... Args>
        requires std::is_base_of_v<Texture, T>
    [[nodiscard]] std::shared_ptr<Texture> GetTexture(const std::string& name, Args&&... args){
        if(name.empty())return nullptr;
        auto it = texture_cache.find(name);
        if(it != texture_cache.end()){
            return it->second;
        }
        return texture_cache.try_emplace(name, std::make_shared<T>(std::forward<Args>(args)...)).first->second;
    }

    template <typename... Args>
    [[nodiscard]] std::shared_ptr<Mesh> getMesh(Args&&... args){
        std::shared_ptr<Mesh> newMesh = std::make_shared<Mesh>(std::forward<Args>(args)...);
        for(auto&& ptr : meshCache){
            if(*ptr == *newMesh){
                return ptr;
            }
        }
        meshCache.push_back(newMesh);
        return newMesh;
    }

    template <typename T, typename... Args>
    std::shared_ptr<Model> CacheModel(const std::string& name, const std::string& path, Args&&... args){
        auto it = modelCache.find(name);
        if(it != modelCache.end()){
            return it->second;
        }
        std::shared_ptr<Model> ptr = std::make_shared<Model>(path);
        ptr->BuildBlas<T>(std::forward<Args>(args)...);
        return modelCache.try_emplace(name, ptr).first->second;
    }

    [[nodiscard]] std::shared_ptr<Model> GetModel(const std::string& name){
        auto it = modelCache.find(name);
        if(it != modelCache.end()){
            return it->second;
        }
        return nullptr;
    }

    void releaseTextures(){
        texture_cache.clear();
    }

    void releaseModels(){
        modelCache.clear();
    }

    void releaseMeshCache(){
        meshCache.clear();
    }
private:
    ResourceManager() = default;
    ~ResourceManager() = default;



    std::vector<std::shared_ptr<Mesh>> meshCache;
    std::unordered_map<std::string, std::shared_ptr<Model>> modelCache;
    std::unordered_map<std::string, std::shared_ptr<Texture>> texture_cache;
};