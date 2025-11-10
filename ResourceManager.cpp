#include "ResourceManager.hpp"
#include "stb_image.h"
#include <iostream>
#include <fstream>


[[nodiscard]] std::shared_ptr<Texture> ResourceManager::GetImageTexture(const std::string& path, bool gammaCorrection, const glm::vec3& colorScale, bool invert){
    if(path.empty())return nullptr;
    auto it = texture_cache.find(path);
    if(it != texture_cache.end()){
        return it->second;
    }
    return texture_cache.try_emplace(path, std::make_shared<ImageTexture>(path, gammaCorrection, colorScale, invert)).first->second;
}


