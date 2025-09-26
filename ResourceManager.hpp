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

    auto load_file(const std::filesystem::path& path) const -> std::string;

    
    auto GetImageTexture(const std::string& path, float tempSRGB = false) -> std::shared_ptr<Texture>;//maybe get image?


    template <typename T, typename... Args>
    requires std::is_base_of_v<Texture,T>
    std::shared_ptr<Texture> GetTexture(const std::string& name,Args&&... args) {
        //instead of giving name, let texture have hash function?
        if(name.empty())return nullptr;
        auto it = texture_cache.find(name);
        if(it != texture_cache.end()){
            return it->second;
        }
        return texture_cache.try_emplace(name,std::make_shared<T>(std::forward<Args>(args)...)).first->second;
    }

    template <typename... Args>
    auto get_model(const std::string& path, Args&&... args) -> std::shared_ptr<Model> {
        auto it = model_cache.find(path);
        if(it == model_cache.end()){
            return model_cache.try_emplace(path,std::make_shared<Model>(path,std::forward<Args>(args)...)).first->second;
        }
        return it->second;
    }
    auto release_textures() -> void;
    private:
    ResourceManager() = default;
    ~ResourceManager() = default;
    std::unordered_map<std::string, std::shared_ptr<Model>> model_cache;
    std::unordered_map<std::string, std::shared_ptr<Texture>> texture_cache;
};