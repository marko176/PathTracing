#include "Texture.hpp"
#include <iostream>

Image::Image(const std::string& filename,float gammaCorrection) {
    data = stbi_load(filename.c_str(),&width,&height,&channels,0);
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

glm::vec3 Image::at(int x,int y) const{
    if(channels == 1) return glm::vec3(data[(y*width + x)],data[(y*width + x)],data[(y*width + x)])/255.0f;
    return glm::vec3{data[(y*width + x)*channels + 0],data[(y*width + x)*channels + 1],data[(y*width + x)*channels + 2]}/255.0f;
}

Image::~Image() {
    stbi_image_free((void*)data);
}


FloatImage::FloatImage(const std::string& filename) {
    data = stbi_loadf(filename.c_str(),&width,&height,&channels,0);
}

glm::vec3 FloatImage::at(int x,int y) const{
    if(channels == 1) return glm::vec3(data[(y*width + x)],data[(y*width + x)],data[(y*width + x)])/255.0f;
    return glm::vec3{data[(y*width + x)*channels + 0],data[(y*width + x)*channels + 1],data[(y*width + x)*channels + 2]}/255.0f;
}



FloatImage::~FloatImage() {
    stbi_image_free((void*)data);
}

SolidColor::SolidColor(const glm::vec3& color) : albedo(color) {}
SolidColor::SolidColor(float r,float g,float b) : albedo(r,g,b) {}

glm::vec3 ImageTexture::texel(int x, int y) const {
    // repeat wrap:
    x = wrap_index(x,image.width);
    y = wrap_index(image.height - y -1,image.height);
    return image.at(x, y);
}



float ImageTexture::alpha(float u,float v) const {
    if(image.channels != 4)return 1;
    float x = u*image.width - 0.5f;
    float y = v*image.height - 0.5f;
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
    /*
    u = glm::fract(u);
    v = 1.0f - glm::fract(v);
    //u = std::clamp(u,0.f,1.f);
    //v = 1.0f - std::clamp(v,0.f,1.f);

    int i = u * image.width;
    int j = std::min<int>(v * image.height,image.height-1);//floor to height-1

    return image.W(i,j);
    */
}

glm::vec3 FloatImageTexture::texel(int x, int y) const {
    // repeat wrap:
    x = wrap_index(x,image.width);
    y = wrap_index(image.height - y -1,image.height);
    return image.at(x, y);
}



float FloatImageTexture::alpha(float u,float v) const {
    if(image.channels != 4)return 1;
    float x = u*image.width - 0.5f;
    float y = v*image.height - 0.5f;
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