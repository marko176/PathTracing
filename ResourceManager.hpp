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

    
    auto get_texture(const std::string& name, float tempSRGB = false) -> std::shared_ptr<Image_texture>;
    auto get_model(const std::string& path) -> std::shared_ptr<Model>;
    auto release_textures() -> void;
    private:
    ResourceManager() = default;
    ~ResourceManager() = default;
    std::unordered_map<std::string, std::shared_ptr<Model>> model_cache;
    std::unordered_map<std::string, std::shared_ptr<Image_texture>> texture_cache;
};