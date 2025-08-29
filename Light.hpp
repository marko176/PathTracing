#pragma once
#include "Hit_record.hpp"
#include "Shape.hpp"
#include <functional>
struct LightSample {
    GeometricInteraction interaction;
    glm::vec3 dir;
};

class Light  {
public:
    virtual bool isDelta() const = 0;
    virtual glm::vec3 L(const GeometricInteraction& interaction, const Ray& ray) const = 0;
    virtual LightSample sample(const glm::vec2& uv) const = 0;
    virtual float PDF(const GeometricInteraction& interaction) const = 0;
    virtual float PDF(const GeometricInteraction& interaction, const Ray& ray) const = 0;
    virtual float Power() const = 0;
    virtual void PreProcess(const AABB& bbox) {}
};

class InfiniteLight : public Light {
public:
    InfiniteLight(const glm::vec3& light_dir, const glm::vec3& light_color,const std::function<float(float)>& powerFunc =  [](float r) -> float { return std::sqrt(r); }) : dir{light_dir} , color{light_color}, powerFunction{powerFunc} {}

    bool isDelta() const final;

    glm::vec3 L(const GeometricInteraction& interaction, const Ray& ray) const override;

    LightSample sample(const glm::vec2& uv) const override;

    float PDF(const GeometricInteraction& interaction) const override;

    float PDF(const GeometricInteraction& interaction, const Ray& ray) const override;

    float Power() const override;

    void PreProcess(const AABB& bbox) override;
private:
    glm::vec3 dir;
    glm::vec3 color;
    float sceneRadius;
    std::function<float(float)> powerFunction;
};

class AreaLight : public Light {
public:
    AreaLight(Shape* light_shape, const glm::vec3& light_color, bool oneSided = false) : shape{light_shape} , color{light_color}, oneSided{oneSided} {}

    bool isDelta() const final;

    glm::vec3 L(const GeometricInteraction& interaction, const Ray& ray) const override ;

    LightSample sample(const glm::vec2& uv) const override ;

    float PDF(const GeometricInteraction& interaction) const override ;

    float PDF(const GeometricInteraction& interaction, const Ray& ray) const override ;

    float Power() const override ;
private:
    Shape* shape;
    glm::vec3 color; //switch to texture/image
    bool oneSided; 
    //add alpha mask
};