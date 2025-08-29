#include "Texture.hpp"
#include <iostream>

Image::Image(const std::string& filename,float tempSRGB) {
    data = stbi_load(filename.c_str(),&width,&height,&channels,0);
    if(tempSRGB && data){
        for(int y = 0;y<height;y++){
            for(int x = 0;x<width;x++){
                for(int i = 0;i<std::min(channels,3);i++){
                    /*
                    unsigned char x = data[(y*width + x)*channels + i];
                    double sRGB = static_cast<int>(x)/255.0;
                    double linear = std::clamp(sRGB_to_linear(sRGB),0.0,1.0);
                    data[(y*width + x)*channels + i] = (unsigned char)std::lround(linear * 255.0);
                    int l = linear * 255.0;
                    std::cout<<(int)sRGBLUT[x]<<" - "<<sRGB<<" "<<(int)x<<" "<<linear<<" "<<l<<"\n";
                    */
                   //watch out about reading uchar as int -> seg fault
                    data[(y*width + x)*channels + i] = sRGBLUT[data[(y*width + x)*channels + i]];
                }
            }
        }
    }
}

glm::vec3 Image::at(int x,int y) const{
    if(data == nullptr)return {1,1,1};
    if(channels == 1) return glm::vec3(data[(y*width + x)],data[(y*width + x)],data[(y*width + x)])/255.0f;
    return glm::vec3{data[(y*width + x)*channels + 0],data[(y*width + x)*channels + 1],data[(y*width + x)*channels + 2]}/255.0f;
}



float Image::W(int x,int y) const{
    if(data == nullptr)return 1;
    if(channels != 4) return 1;
    return data[(y*width + x)*4 + 3]/255.0f;
}

Image::~Image() {
    if(data != nullptr)
    stbi_image_free((void*)data);
}





Solid_color::Solid_color(const glm::vec3& color) : albedo(color) {}
Solid_color::Solid_color(float r,float g,float b) : albedo(r,g,b) {}
glm::vec3 Solid_color::color_value(float u,float v) const {
    return albedo;
}









Image_texture::Image_texture(const std::string& filename,float tempSRGB) : image(filename,tempSRGB) {}

glm::vec3 Image_texture::color_value(float u, float v) const {
    //glm::fract?
    u = glm::fract(u);
    v = 1.0f - glm::fract(v);
    //u = std::clamp(u,0.f,1.f);
    //v = 1.0f - std::clamp(v,0.f,1.f);

    int i = u * image.width;
    int j = std::min<int>(v * image.height,image.height-1);//floor to height-1

    return image.at(i,j);
}

inline int wrap_index(int i, int n) {
    int m = i % n;
    if (m < 0) m += n;
    return m;
}

glm::vec3 Image_texture::texel(int x, int y) const {
    if (!image.data) return glm::vec3(1.0f);
    // repeat wrap:
    x = wrap_index(x,image.width);
    y = wrap_index(image.height - y -1,image.height);
    return image.at(x, y);
}



float Image_texture::alpha(float u,float v) const {
    u = glm::fract(u);
    v = 1.0f - glm::fract(v);
    //u = std::clamp(u,0.f,1.f);
    //v = 1.0f - std::clamp(v,0.f,1.f);

    int i = u * image.width;
    int j = std::min<int>(v * image.height,image.height-1);//floor to height-1

    return image.W(i,j);
}
