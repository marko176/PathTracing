#pragma once
#include "Interaction.hpp"
#include "Shape.hpp"
#include <functional>

struct LightSample {
    GeometricInteraction interaction;
    glm::vec3 dir;
};

class Light  {
public:
    virtual ~Light() = default;
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
    virtual ~InfiniteLight() = default;

    InfiniteLight(const glm::vec3& light_dir, const glm::vec3& light_color,const std::function<float(float)>& powerFunc = [](float r) -> float { return std::sqrt(r); }) : dir{light_dir} , color{light_color}, powerFunction{powerFunc} {}

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

class PointLight : public Light {
public:
    virtual ~PointLight() = default;

    PointLight(const glm::vec3& p, const glm::vec3& light_color,const std::function<float(float)>& powerFunc = [](float r) -> float { return 4 * r; }) : p{p} , color{light_color}, powerFunction{powerFunc} {}

    bool isDelta() const final;

    glm::vec3 L(const GeometricInteraction& interaction, const Ray& ray) const override;

    LightSample sample(const glm::vec2& uv) const override;

    float PDF(const GeometricInteraction& interaction) const override;

    float PDF(const GeometricInteraction& interaction, const Ray& ray) const override;

    float Power() const override;

    void PreProcess(const AABB& bbox) override;
private:
    glm::vec3 p;
    glm::vec3 color;
    float sceneRadius;
    std::function<float(float)> powerFunction;
};

class AreaLight : public Light {
public:
    virtual ~AreaLight() = default;

    AreaLight(const std::shared_ptr<Shape>& light_shape, const glm::vec3& light_color, bool oneSided = false) : shape{light_shape} , color{light_color}, oneSided{oneSided} {}

    bool isDelta() const final;

    glm::vec3 L(const GeometricInteraction& interaction, const Ray& ray) const override ;

    LightSample sample(const glm::vec2& uv) const override ;

    float PDF(const GeometricInteraction& interaction) const override ;

    float PDF(const GeometricInteraction& interaction, const Ray& ray) const override ;

    float Power() const override ;

    //tranform -> changes shape to tranformed shape 
private:
    std::shared_ptr<Shape> shape;
    glm::vec3 color; //switch to texture/image
    bool oneSided; 
    //add alpha mask
};

class TransformedLight : public Light {
public:
    virtual ~TransformedLight() = default;

    TransformedLight(const glm::mat4& transform,const std::shared_ptr<Light>& light) : transform(transform), light(light), normalMatrix(glm::transpose(glm::inverse(glm::mat3(transform)))), invTransform(glm::inverse(transform)) {}

    bool isDelta() const final;

    glm::vec3 L(const GeometricInteraction& interaction, const Ray& ray) const override ;

    LightSample sample(const glm::vec2& uv) const override ;

    float PDF(const GeometricInteraction& interaction) const override ;

    float PDF(const GeometricInteraction& interaction, const Ray& ray) const override ;

    float Power() const override ;

private:
    glm::mat4 transform;
    std::shared_ptr<Light> light;
    glm::mat3 normalMatrix;
    glm::mat4 invTransform;
};