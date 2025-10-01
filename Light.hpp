#pragma once
#include "Interaction.hpp"
#include "Shape.hpp"
#include <functional>

struct LightSample {
    //should have light color!
    glm::vec3 L;
    GeometricInteraction interaction;
    glm::vec3 dir;
};

class Light  {
public:
    virtual ~Light() = default;
    virtual bool isDelta() const = 0;
    virtual glm::vec3 L(const GeometricInteraction& interaction, const Ray& ray) const = 0;
    virtual LightSample sample(const glm::vec2& uv, float time) const = 0;
    virtual float PDF(const GeometricInteraction& interaction, float time) const = 0;
    virtual float PDF(const GeometricInteraction& interaction, const Ray& ray) const = 0;
    virtual float Power() const = 0;
    virtual void PreProcess(const AABB& bbox) {}
};

//infinite light has Le which doesnt need interaction!
//faraway / infinite directional
//uniform infinite
//texture infinite

class InfiniteLight : public Light {
public:
    virtual ~InfiniteLight() = default;

    virtual glm::vec3 Le(const Ray& ray) const = 0;
    bool isDelta() const override;
    float PDF(const GeometricInteraction& interaction, float time) const override;
    float PDF(const GeometricInteraction& interaction, const Ray& ray) const override;
    void PreProcess(const AABB& bbox) override;
protected:
    float sceneRadius;
};

class UniformInfiniteLight : public InfiniteLight {
public:
    virtual ~UniformInfiniteLight() = default;

    UniformInfiniteLight(const glm::vec3& light_color,const std::function<float(float)>& powerFunc = [](float r) -> float { return std::sqrt(r); }) : color{light_color}, powerFunction{powerFunc} {}

    glm::vec3 Le(const Ray& ray) const override;

    glm::vec3 L(const GeometricInteraction& interaction, const Ray& ray) const override;

    LightSample sample(const glm::vec2& uv, float time) const override;

    float Power() const override;

private:
    glm::vec3 color;
    std::function<float(float)> powerFunction;
};

class FunctionInfiniteLight : public InfiniteLight {
public:
    virtual ~FunctionInfiniteLight() = default;

    FunctionInfiniteLight(const std::function<glm::vec3(const Ray& ray)>& lightFunc,const std::function<float(float)>& powerFunc = [](float r) -> float { return std::sqrt(r); }) : lightFunction(lightFunc), powerFunction{powerFunc} {}

    glm::vec3 Le(const Ray& ray) const override;

    glm::vec3 L(const GeometricInteraction& interaction, const Ray& ray) const override;

    LightSample sample(const glm::vec2& uv, float time) const override;

    float Power() const override;

private:
    std::function<glm::vec3(const Ray& ray)> lightFunction;
    std::function<float(float)> powerFunction;
};

class TextureInfiniteLight : public InfiniteLight {
public:
    virtual ~TextureInfiniteLight() = default;

    TextureInfiniteLight(const std::shared_ptr<Texture>& tex, float LeScale,const std::function<float(float)>& powerFunc = [](float r) -> float { return std::sqrt(r); }) : tex{tex}, LeScale(LeScale), powerFunction{powerFunc} {}

    glm::vec3 Le(const Ray& ray) const override;

    glm::vec3 L(const GeometricInteraction& interaction, const Ray& ray) const override;

    LightSample sample(const glm::vec2& uv, float time) const override;

    float Power() const override;

private:
    std::shared_ptr<Texture> tex;
    float LeScale;
    std::function<float(float)> powerFunction;
};


//distantLight : InfiniteLIght?
class DistantLight : public Light {
public:
    virtual ~DistantLight() = default;

    DistantLight(const glm::vec3& light_dir, const glm::vec3& light_color,const std::function<float(float)>& powerFunc = [](float r) -> float { return std::sqrt(r); }) : dir{light_dir} , color{light_color}, powerFunction{powerFunc} {}

    bool isDelta() const final;

    glm::vec3 L(const GeometricInteraction& interaction, const Ray& ray) const override;

    LightSample sample(const glm::vec2& uv, float time) const override;

    float PDF(const GeometricInteraction& interaction, float time) const override;

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

    LightSample sample(const glm::vec2& uv, float time) const override;

    float PDF(const GeometricInteraction& interaction, float time) const override;

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

    LightSample sample(const glm::vec2& uv, float time) const override ;

    float PDF(const GeometricInteraction& interaction, float time) const override ;

    float PDF(const GeometricInteraction& interaction, const Ray& ray) const override ;

    float Power() const override ;

    std::shared_ptr<Shape> getShape() const ;
private:
    std::shared_ptr<Shape> shape;
    glm::vec3 color; //switch to texture/image
    bool oneSided; 
    //add alpha mask
};

class TransformedLight : public Light {
public:
    virtual ~TransformedLight() = default;

    TransformedLight(const std::shared_ptr<Light>& light,const glm::mat4& transform) : light(light), transform(transform), normalMatrix(glm::transpose(glm::inverse(glm::mat3(transform)))), invTransform(glm::inverse(transform)) {}

    bool isDelta() const final;

    glm::vec3 L(const GeometricInteraction& interaction, const Ray& ray) const override ;

    LightSample sample(const glm::vec2& uv, float time) const override ;

    float PDF(const GeometricInteraction& interaction, float time) const override ;

    float PDF(const GeometricInteraction& interaction, const Ray& ray) const override ;

    float Power() const override ;

private:
    std::shared_ptr<Light> light;
    glm::mat4 transform;
    glm::mat3 normalMatrix;
    glm::mat4 invTransform;
};

class AnimatedLight : public Light {
public:
    virtual ~AnimatedLight() = default;

    AnimatedLight(const std::shared_ptr<Light>& light,const glm::vec3& direction, const glm::vec2& timeBounds) : light(light), dir(direction), timeBounds(timeBounds) {}

    bool isDelta() const final;

    glm::vec3 L(const GeometricInteraction& interaction, const Ray& ray) const override ;

    LightSample sample(const glm::vec2& uv, float time) const override ;

    float PDF(const GeometricInteraction& interaction, float time) const override ;

    float PDF(const GeometricInteraction& interaction, const Ray& ray) const override ;

    float Power() const override ;

private:
    std::shared_ptr<Light> light;
    glm::vec3 dir;
    glm::vec2 timeBounds;
};