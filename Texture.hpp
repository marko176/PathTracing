#pragma once
#include <memory>
#include "stb_image.h"
#include <iostream>
#include <array>
#include "Interaction.hpp"

inline int wrap_index(int i, int n){
    int m = i % n;
    if(m < 0) m += n;
    return m;
}

template <std::floating_point T>
inline constexpr T linear_to_sRGB(T linear){
    linear = glm::clamp<T>(linear, 0.0, 1.0);
    return linear < static_cast<T>(0.0031308) ? static_cast<T>(12.92) * linear : static_cast<T>(1.055) * std::pow(linear, 1.0 / 2.4) - 0.055;
}

template <std::floating_point T>
inline constexpr T sRGB_to_linear(T sRGB){
    sRGB = glm::clamp<T>(sRGB, 0.0, 1.0);
    if(sRGB <= static_cast<T>(0.04045)) return sRGB / static_cast<T>(12.92);
    return std::pow<T>((sRGB + 0.055) / 1.055, 2.4);
}

inline const std::array<unsigned char, 256> sRGBLUT = [](){
    std::array<unsigned char, 256> tempLUT;
    for(int i = 0;i < 256;i++){
        double sRGB = static_cast<unsigned char>(i) / 255.0;
        double linear = glm::clamp(sRGB_to_linear(sRGB), 0.0, 1.0);
        tempLUT[i] = static_cast<unsigned char>(std::lround(linear * 255.0));
    }
    return tempLUT;
    }();


//add FloatImage -> for alpha mask -> no division needed

struct Image{

    Image(const std::string_view filename, bool gammaCorrection = false);
    //getChannel
    float GetChannelAt(const glm::ivec2& p, int ch) const{//pass wrap mode
        //channel 1,2,3,4
        int x = wrap_index(p.x, width);
        int y = wrap_index(p.y, height);
        return data[(y * width + x) * channels + (ch - 1)] / 255.0f;
    }

    float NearestChannel(const glm::vec2& p, int ch) const{
        glm::ivec2 pi = { p.x * width, p.y * height };
        return GetChannelAt(pi, ch);
    }

    glm::ivec2 Resolution() const{
        return { width,height };
    }

    int Channels() const{
        return channels;
    }



    ~Image();
private:
    unsigned char* data;
    int width;
    int height;
    int channels;
};

struct FloatImage{


    FloatImage(const std::string_view filename);
    //getChannel
    float GetChannelAt(const glm::ivec2& p, int ch) const{//pass wrap mode
        //channel 1,2,3,4
        int x = wrap_index(p.x, width);
        int y = wrap_index(p.y, height);
        return data[(y * width + x) * channels + (ch - 1)];
    }

    float NearestChannel(const glm::vec2& p, int ch) const{
        glm::ivec2 pi = { p.x * width, p.y * height };
        return GetChannelAt(pi, ch);
    }

    glm::ivec2 Resolution() const{
        return { width,height };
    }

    int Channels() const{
        return channels;
    }



    ~FloatImage();
private:
    float* data;
    int width;
    int height;
    int channels;
};

class Texture{
public:
    Texture(const glm::vec3& colorScale = glm::vec3(1, 1, 1), bool invert = false) : colorScale(colorScale), invert(invert){}
    virtual ~Texture() = default;
    virtual float alpha(float u, float v) const{
        return 1;
    }
    virtual int Channels() const = 0;
    virtual glm::vec3 Evaluate(const SurfaceInteraction& interaction) const = 0;
protected:
    glm::vec3 colorScale;
    bool invert;
};

class SolidColor : public Texture{
public:
    virtual ~SolidColor() = default;
    SolidColor(const glm::vec3& color, const glm::vec3& colorScale = glm::vec3(1), bool invert = false);
    SolidColor(float r, float g, float b, const glm::vec3& colorScale = glm::vec3(1), bool invert = false);

    glm::vec3 Evaluate(const SurfaceInteraction& interaction) const override{
        return colorScale * albedo;
    }
    virtual int Channels() const{ return 3; }
private:
    glm::vec3 albedo;
};

class ImageTexture : public Texture{
public:
    virtual ~ImageTexture() = default;
    ImageTexture(const std::string& filename, bool gammaCorrection = false, const glm::vec3& colorScale = glm::vec3(1), bool invert = false) : Texture(colorScale, invert), image(filename, gammaCorrection){};

    float alpha(float u, float v) const override;

    glm::vec3 Evaluate(const SurfaceInteraction& interaction) const override{ //take enum? bilerp
        glm::ivec2 res = image.Resolution();
        float x = interaction.uv.x * res.x - 0.5f;
        float y = interaction.uv.y * res.y - 0.5f;
        int xi = std::floor(x);
        int yi = std::floor(y);
        float dx = x - xi;
        float dy = y - yi;

        glm::vec3 a = texel(xi, yi);//image.texel
        glm::vec3 b = texel(xi + 1, yi);
        glm::vec3 c = texel(xi, yi + 1);
        glm::vec3 d = texel(xi + 1, yi + 1);
        return colorScale * ((1 - dx) * (1 - dy) * a + dx * (1 - dy) * b +
            (1 - dx) * dy * c + dx * dy * d);
    }
    virtual int Channels() const{
        return image.Channels();
    }
private:
    glm::vec3 texel(int x, int y) const;
    Image image;
};

class FloatImageTexture : public Texture{
public:
    virtual ~FloatImageTexture() = default;
    FloatImageTexture(const std::string& filename, const glm::vec3& colorScale = glm::vec3(1), bool invert = false) : Texture(colorScale, invert), image(filename){};

    float alpha(float u, float v) const override;

    glm::vec3 Evaluate(const SurfaceInteraction& interaction) const override{
        glm::ivec2 res = image.Resolution();
        float x = interaction.uv.x * res.x - 0.5f;
        float y = interaction.uv.y * res.y - 0.5f;
        int xi = std::floor(x);
        int yi = std::floor(y);
        float dx = x - xi;
        float dy = y - yi;

        glm::vec3 a = texel(xi, yi);//image.texel
        glm::vec3 b = texel(xi + 1, yi);
        glm::vec3 c = texel(xi, yi + 1);
        glm::vec3 d = texel(xi + 1, yi + 1);
        return colorScale * ((1 - dx) * (1 - dy) * a + dx * (1 - dy) * b +
            (1 - dx) * dy * c + dx * dy * d);
    }
    virtual int Channels() const{ return image.Channels(); }
private:
    glm::vec3 texel(int x, int y) const;
    FloatImage image;
};

class CheckerTexture : public Texture{
public:
    virtual ~CheckerTexture() = default;
    CheckerTexture(const std::shared_ptr<Texture>& textureA, const std::shared_ptr<Texture>& textureB, const glm::vec2& uvscale, const glm::vec3& colorScale = glm::vec3(1), bool invert = false) : Texture(colorScale, invert), tex1(textureA), tex2(textureB), invScale(1.0f / uvscale){}
    float alpha(float u, float v) const override;


    glm::vec3 Evaluate(const SurfaceInteraction& interaction) const override{
        glm::ivec2 uv = glm::floor(interaction.uv * invScale);
        if((uv.x + uv.y) % 2 == 0)return colorScale * tex1->Evaluate(interaction);
        return colorScale * tex2->Evaluate(interaction);
    }
    virtual int Channels() const{ return tex1->Channels(); }
private:
    std::shared_ptr<Texture> tex1;
    std::shared_ptr<Texture> tex2;
    glm::vec2 invScale;
};

class UVTexture : public Texture{
public:
    virtual ~UVTexture() = default;
    glm::vec3 Evaluate(const SurfaceInteraction& interaction) const override{
        return colorScale * glm::vec3 { interaction.uv,0 };
    }
    virtual int Channels() const{ return 3; }
};

class NormalTexture : public Texture{
public:
    virtual ~NormalTexture() = default;
    glm::vec3 Evaluate(const SurfaceInteraction& interaction) const override{
        return colorScale * interaction.ns;
    }
    virtual int Channels() const{ return 3; }
};





/*
template <typename T>
class Texture2 {
public:
    virtual T Evaluate(const SurfaceInteraction& interaction) const = 0;
    virtual ~Texture2() = default;
};

template <typename T>
class SolidTexture : public Texture2<T>{
public:
    SolidTexture(const T& val) : value{val} {}
    T Evaluate(const SurfaceInteraction& interaction) const{
        return value;
    }
protected:
    T value;
};

template <typename T, typename U>
class ConvolutedTexture : public Texture2<U>{
public:
    ConvolutedTexture(const std::shared_ptr<Texture2<T>>& tex1, const std::shared_ptr<Texture2<U>>& tex2) : texture1{tex1}, texture2{tex2} {}
    U Evaluate(const SurfaceInteraction& interaction) const{
        return texture1->Evaluate(interaction) * texture2->Evaluate(interaction);
    }
protected:
    std::shared_ptr<Texture<T>> texture1;
    std::shared_ptr<Texture<U>> texture2;
};

template <typename T>
class ImageTexture2 : public Texture2<T>{
public:
    ImageTexture(const std::string& file, bool gammaCorrect = false, bool invert = false, int activeChannelMask = 0b1111) {
        //mask says what channels to take

        if(invert)
            stbi_set_flip_vertically_on_load(true);

        unsigned char* tempData = stbi_load(file.data(),&width,&height,&channels,4);
        if(tempData == nullptr){
            std::cerr<<"Failed to load image: "<<file<<std::endl;
            std::abort();
        }
        if(gammaCorrect){
            for(int y = 0;y<height;y++){
                for(int x = 0;x<width;x++){
                    for(int i = 0;i<3;i++){//3 becouse never do gammaCorrection on alpha
                        tempData[(y*width + x)*channels + i] = sRGBLUT[tempData[(y*width + x)*channels + i]];
                    }
                }
            }
        }
        int activeChannels =    ((activeChannelMask & 0b1000) >> 3) +
                                ((activeChannelMask & 0b0100) >> 2) +
                                ((activeChannelMask & 0b0010) >> 1) +
                                ((activeChannelMask & 0b0001) >> 0);
        data = new unsigned char[width*height*activeChannels];
        for(int y = 0;y<height;y++){
            for(int x = 0;x<width;x++){
                for(int i = 0;i<activeChannels;i++){//3 becouse never do gammaCorrection on alpha
                    data[(y*width + x)*activeChannels + i];
                }
            }
        }
        channels = activeChannels;
        if(invert)
            stbi_set_flip_vertically_on_load(false);
        stbi_image_free((void*)tempData);
    }

    ~ImageTexture() {
        stbi_image_free((void*)data);
    }

    T Evaluate(const SurfaceInteraction& interaction) const{
        return value;
    }
protected:

    unsigned char* data;
    int width;
    int height;
    int channels;
};

template <typename T>
class HDRImageTexture : public Texture2<T>{
public:
    ~HDRImageTexture() {
        stbi_image_free((void*)data);
    }
protected:
    float* data;
    int width;
    int height;
    int channels;
};


*/
