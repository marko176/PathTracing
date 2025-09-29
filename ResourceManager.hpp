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
    
    auto GetImageTexture(const std::string& path, float tempSRGB = false) -> std::shared_ptr<Texture>;//maybe get image?


    template <typename T, typename... Args>
    requires std::is_base_of_v<Texture,T>
    std::shared_ptr<Texture> GetTexture(const std::string& name,Args&&... args) {
        if(name.empty())return nullptr;
        auto it = texture_cache.find(name);
        if(it != texture_cache.end()){
            return it->second;
        }
        return texture_cache.try_emplace(name,std::make_shared<T>(std::forward<Args>(args)...)).first->second;
    }

    template <typename... Args>
    auto getMesh(Args&&... args) -> std::shared_ptr<Mesh> {
        std::shared_ptr<Mesh> newMesh = std::make_shared<Mesh>(std::forward<Args>(args)...);
        for(auto&& ptr : meshCache){
            if(*ptr == *newMesh){
                return ptr;
            }
        }
        meshCache.push_back(newMesh);
        return newMesh;
    }

    template <typename... Args>
    auto GetModel(const std::string& name, Args&&... args) -> std::shared_ptr<Model> {
        auto it = modelCache.find(name);
        if(it != modelCache.end()){
            return it->second;
        }
        return modelCache.try_emplace(name,std::make_shared<Model>(std::forward<Args>(args)...)).first->second;
    }

    void releaseTextures();
private:
    ResourceManager() = default;
    ~ResourceManager() = default;



    std::vector<std::shared_ptr<Mesh>> meshCache;
    std::unordered_map<std::string, std::shared_ptr<Model>> modelCache;
    std::unordered_map<std::string, std::shared_ptr<Texture>> texture_cache;
};