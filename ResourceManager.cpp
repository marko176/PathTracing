#include "ResourceManager.hpp"
#include "stb_image.h"
#include <iostream>
#include <fstream>

auto ResourceManager::load_file(const std::filesystem::path& path) const -> std::string {
    std::ifstream in(path, std::ios::in);
    return std::string(std::istreambuf_iterator<char>(in),std::istreambuf_iterator<char>());
}



auto ResourceManager::get_texture(const std::string& path, float tempSRGB) -> std::shared_ptr<ImageTexture> {
    if(path.empty())return nullptr;
    auto it = texture_cache.find(path);
    if(it != texture_cache.end()){
        return it->second;
    }
    return texture_cache.try_emplace(path,std::make_shared<ImageTexture>(path,tempSRGB)).first->second;
}



auto ResourceManager::release_textures() -> void {
    
}