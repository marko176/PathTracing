#include "ResourceManager.hpp"
#include "stb_image.h"
#include <iostream>
#include <fstream>


auto ResourceManager::GetImageTexture(const std::string& path, float tempSRGB) -> std::shared_ptr<Texture> {
    if(path.empty())return nullptr;
    auto it = texture_cache.find(path);
    if(it != texture_cache.end()){
        return it->second;
    }
    return texture_cache.try_emplace(path,std::make_shared<ImageTexture>(path,tempSRGB)).first->second;
}



void ResourceManager::releaseTextures() {
    texture_cache.clear();
}