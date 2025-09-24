#pragma once
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtx/intersect.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <memory>
#include "stb_image.h"
#include <iostream>
#include <array>
#include "Interaction.hpp"


inline constexpr double linear_to_sRGB(double linear){
    if(linear != linear)std::cout<<"linera_component is NaN\n";
    linear = std::clamp(linear,0.0,1.0);
    if(linear > 0){
        return linear < 0.0031308 ? 12.92 * linear : 1.055 * std::pow(linear,1.0/2.4) - 0.055;
    }
    return 0;
}

inline constexpr double sRGB_to_linear(double sRGB){
    sRGB = std::clamp(sRGB,0.0,1.0);
    if (sRGB <= 0.04045) return sRGB / 12.92;
    return std::pow((sRGB + 0.055) / 1.055, 2.4);
}

static inline std::array<unsigned char,256> sRGBLUT = [](){
    std::array<unsigned char,256> LUT;
    for(int i = 0;i<256;i++){
        double sRGB = static_cast<unsigned char>(i)/255.0;
        double linear = std::clamp(sRGB_to_linear(sRGB),0.0,1.0);
        LUT[i] = (unsigned char)std::lround(linear * 255.0);
    }
    return LUT;
}();


//add FloatImage -> for alpha mask -> no division needed

struct Image{
    unsigned char* data;
    int width;
    int height;
    int channels;

    Image(const std::string& filename,float gammaCorrection = false);
    glm::vec3 at(int x,int y) const;
    //getChannel
    float GetChannelAt(const glm::ivec2& p,int ch) const {
        return data[(p.y*width + p.x)*channels + ch]/255.0f;
    }
    float W(int x,int y) const;
    ~Image() ;

};

class Texture{
public:
    virtual ~Texture() = default;
    virtual float alpha(float u,float v) const {
        return 1;
    }

    virtual glm::vec3 Evaluate(const SurfaceInteraction& interaction) const = 0;
};

class SolidColor : public Texture {
public:
    virtual ~SolidColor() = default;
    SolidColor(const glm::vec3& color);
    SolidColor(float r,float g,float b);


    glm::vec3 Evaluate(const SurfaceInteraction& interaction) const override {
        return albedo;
    }
private:
    glm::vec3 albedo;
};

class ImageTexture : public Texture {
public:
    virtual ~ImageTexture() = default;
    ImageTexture(const std::string& filename,float gammaCorrection = false);

    float alpha(float u,float v) const override;
    
    glm::vec3 Evaluate(const SurfaceInteraction& interaction) const override { //take enum? bilerp
        float x = interaction.uv.x*image.width - 0.5f;
        float y = interaction.uv.y*image.height - 0.5f;
        int xi = std::floor(x);
        int yi = std::floor(y);
        float dx = x - xi;
        float dy = y - yi;
        
        glm::vec3 a = texel(xi,yi);//image.texel
        glm::vec3 b = texel(xi+1,yi);
        glm::vec3 c = texel(xi,yi+1);
        glm::vec3 d = texel(xi+1,yi+1);
        return ((1 - dx) * (1 - dy) * a + dx * (1 - dy) * b +
        (1 - dx) *      dy  * c + dx *      dy  * d);
    }
    
private: 
    glm::vec3 texel(int x, int y) const;
    Image image;
};

class CheckerTexture : public Texture {
public:
    virtual ~CheckerTexture() = default;
    CheckerTexture(const std::shared_ptr<Texture>& tex1, const std::shared_ptr<Texture>& tex2,const glm::vec2& scale) : tex1(tex1), tex2(tex2), invScale(1.0f/scale) {}
    float alpha(float u,float v) const override;
    

    glm::vec3 Evaluate(const SurfaceInteraction& interaction) const override{
        glm::ivec2 uv = glm::floor(interaction.uv * invScale);
        if((uv.x+uv.y)%2 == 0)return tex1->Evaluate(interaction);
        return tex2->Evaluate(interaction);
    }
private:
    std::shared_ptr<Texture> tex1;
    std::shared_ptr<Texture> tex2;
    glm::vec2 invScale;
};

class UVTexture : public Texture {
public:
    virtual ~UVTexture() = default;
    glm::vec3 Evaluate(const SurfaceInteraction& interaction) const override{
        return {interaction.uv,0};
    }
};

class NormalTexture : public Texture {
public:
    virtual ~NormalTexture() = default;
    glm::vec3 Evaluate(const SurfaceInteraction& interaction) const override{
        return interaction.ns;
    }
};