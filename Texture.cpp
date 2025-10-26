#include "Texture.hpp"
#include <iostream>

Image::Image(const std::string_view filename,bool gammaCorrection) {
    data = stbi_load(filename.data(),&width,&height,&channels,0);
    if(data == nullptr){
        std::cerr<<"Failed to load image: "<<filename<<std::endl;
        std::abort();
    }
    if(gammaCorrection && data){
        for(int y = 0;y<height;y++){
            for(int x = 0;x<width;x++){
                for(int i = 0;i<std::min(channels,3);i++){
                    data[(y*width + x)*channels + i] = sRGBLUT[data[(y*width + x)*channels + i]];
                }
            }
        }
    }
}

Image::~Image() {
    stbi_image_free((void*)data);
}


FloatImage::FloatImage(const std::string_view filename) {
    data = stbi_loadf(filename.data(),&width,&height,&channels,0);
    if(data == nullptr){
        std::cerr<<"Failed to load image: "<<filename<<std::endl;
        std::abort();
    }
}

FloatImage::~FloatImage() {
    stbi_image_free((void*)data);
}

SolidColor::SolidColor(const glm::vec3& color) : albedo(color) {}
SolidColor::SolidColor(float r,float g,float b) : albedo(r,g,b) {}

glm::vec3 ImageTexture::texel(int x, int y) const {
    glm::ivec2 p = {x,y};
    return {image.GetChannelAt(p,1),image.GetChannelAt(p,2),image.GetChannelAt(p,3)};
}



float ImageTexture::alpha(float u,float v) const {
    if(image.Channels() != 4)return 1;
    glm::ivec2 res = image.Resolution();
    float x = u*res.x - 0.5f;
    float y = v*res.y - 0.5f;
    int xi = std::floor(x);
    int yi = std::floor(y);
    float dx = x - xi;
    float dy = y - yi;
    
    float a = image.GetChannelAt({xi,yi},4);//image.image.GetChannelAt
    float b = image.GetChannelAt({xi+1,yi},4);
    float c = image.GetChannelAt({xi,yi+1},4);
    float d = image.GetChannelAt({xi+1,yi+1},4);
    return ((1 - dx) * (1 - dy) * a + dx * (1 - dy) * b +
    (1 - dx) *      dy  * c + dx *      dy  * d);
}

glm::vec3 FloatImageTexture::texel(int x, int y) const {
    // repeat wrap:
    glm::ivec2 p = {x,y};
    return {image.GetChannelAt(p,1),image.GetChannelAt(p,2),image.GetChannelAt(p,3)};
}



float FloatImageTexture::alpha(float u,float v) const {
    if(image.Channels() != 4)return 1;
    glm::ivec2 res = image.Resolution();
    float x = u*res.x - 0.5f;
    float y = v*res.y - 0.5f;
    int xi = std::floor(x);
    int yi = std::floor(y);
    float dx = x - xi;
    float dy = y - yi;
    
    float a = image.GetChannelAt({xi,yi},4);//image.image.GetChannelAt
    float b = image.GetChannelAt({xi+1,yi},4);
    float c = image.GetChannelAt({xi,yi+1},4);
    float d = image.GetChannelAt({xi+1,yi+1},4);
    return ((1 - dx) * (1 - dy) * a + dx * (1 - dy) * b +
    (1 - dx) *      dy  * c + dx *      dy  * d);
}

float CheckerTexture::alpha(float u,float v) const {
    glm::ivec2 uv = glm::floor(glm::vec2{u,v} * invScale);
    if((uv.x+uv.y)%2 == 0)return tex1->alpha(u,v);
    return tex2->alpha(u,v);
}