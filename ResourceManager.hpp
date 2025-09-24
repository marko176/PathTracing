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

    
    auto get_texture(const std::string& name, float tempSRGB = false) -> std::shared_ptr<ImageTexture>;
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
    std::unordered_map<std::string, std::shared_ptr<ImageTexture>> texture_cache;
};