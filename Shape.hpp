#pragma once
#include <memory>
#include "Ray.hpp"
#include "AABB.hpp"
#include "Material.hpp"
#include <mutex>
#include <numbers>
class Mesh;
class Shape{
public:
    virtual AABB BoundingBox() const = 0;
    virtual bool IntersectPred(const Ray& ray, float max) const = 0;
    virtual bool Intersect(const Ray& ray, SurfaceInteraction& interaction, float max) const = 0;
    virtual float Area() const = 0;
    virtual SurfaceInteraction Sample(const glm::vec2& u) const = 0;
    virtual float PDF(const GeometricInteraction& interaction) const = 0;
    virtual float PDF(const GeometricInteraction& interaction, const Ray& ray) const = 0; //should just give normal info (shape sample context, maybe get that from sample?) 
};

class SphereShape : public Shape{
public:
    SphereShape(const glm::vec3& sphereCenter, float sphereRadius) : center(sphereCenter), bbox {}, radius(sphereRadius){
        glm::vec3 rvec = glm::vec3(radius);
        bbox.Expand(center - rvec);
        bbox.Expand(center + rvec);
    }
    bool Intersect(const Ray& ray, SurfaceInteraction& interaction, float max) const override;

    bool IntersectPred(const Ray& ray, float max) const override;

    AABB BoundingBox() const override{
        return bbox;
    }

    static glm::vec2 getSphereUV(glm::vec3 p){
        p = glm::normalize(p);
        float theta = std::acos(glm::clamp(p.y, -1.0f, 1.0f));
        float phi = std::atan2(p.z, p.x);
        if(phi < 0) phi += 2.0f * std::numbers::pi_v<float>;
        float u = std::numbers::inv_pi_v<float> *phi * 0.5f;
        float v = std::numbers::inv_pi_v<float> *theta;
        return { u,v };
    }


    SurfaceInteraction Sample(const glm::vec2& u) const override;

    float Area() const override;

    float PDF(const GeometricInteraction& interaction) const override;

    float PDF(const GeometricInteraction& interaction, const Ray& ray) const override;

private:
    glm::vec3 center;
    AABB bbox;
    float radius;
};

class TriangleShape : public Shape{
public:
    TriangleShape(uint32_t meshIndex, uint32_t triangleIndex) : MeshIndex(meshIndex), TriIndex(triangleIndex){}

    bool Intersect(const Ray& ray, SurfaceInteraction& interaction, float max) const override;

    bool IntersectPred(const Ray& ray, float max) const override;

    AABB BoundingBox() const override;

    SurfaceInteraction Sample(const glm::vec2& u) const override;

    float Area() const override;

    float PDF(const GeometricInteraction& interaction) const override;

    float PDF(const GeometricInteraction& interaction, const Ray& ray) const override;

    //maybe move this to resource manager?
    static uint32_t addMesh(Mesh* mesh){
        const std::lock_guard<std::mutex> ml(meshListLock);
        for(std::size_t i = 0;i < meshList.size();i++){
            if(meshList[i] == nullptr || meshList[i] == mesh){
                meshList[i] = mesh;
                return i;
            }
        }
        meshList.push_back(mesh);
        return static_cast<uint32_t>(meshList.size() - 1);
    }

    static void removeMesh(Mesh* mesh){
        const std::lock_guard<std::mutex> ml(meshListLock);
        for(std::size_t i = 0;i < meshList.size();i++){
            if(meshList[i] == mesh){
                meshList[i] = nullptr;
            }
        }
    }
private:
    uint32_t MeshIndex;
    uint32_t TriIndex;
    static inline std::vector<Mesh*> meshList;
    static inline std::mutex meshListLock;
};

class QuadShape : public Shape{
public:
    QuadShape(const glm::vec3& origin, const glm::vec3& u, const glm::vec3& v) : Q(origin), u(u), v(v), bbox {}{
        glm::vec3 n = glm::cross(u, v);
        normal = glm::normalize(n);
        D = glm::dot(normal, Q);
        bbox.Expand(Q);
        bbox.Expand(Q + u + v);
        bbox.Expand(Q + u);
        bbox.Expand(Q + v);
        w = n / glm::dot(n, n);
    }

    bool Intersect(const Ray& ray, SurfaceInteraction& interaction, float max) const override;

    bool IntersectPred(const Ray& ray, float max) const override;

    AABB BoundingBox() const override{
        return bbox;
    }

    SurfaceInteraction Sample(const glm::vec2& uv) const override{
        return GeometricInteraction { Q + uv.x * u + uv.y * v, normal };
    }

    float Area() const override{
        return glm::length(glm::cross(u, v));
    }

    float PDF(const GeometricInteraction& interaction) const override{
        return 1.0f / Area();
    }

    float PDF(const GeometricInteraction& interaction, const Ray& ray) const override{
        glm::vec3 to_shape = interaction.p - ray.origin;
        float dist_squared = glm::dot(to_shape, to_shape);
        float light_cosine = std::abs(glm::dot(-ray.dir, interaction.n));
        float area = Area();
        if(area == 0)return 0;
        return dist_squared / (light_cosine * area);
    }
private:

    static bool is_interior(float a, float b){
        return a >= 0 && a <= 1 && b >= 0 && b <= 1;
    }
    glm::vec3 Q;
    glm::vec3 u;
    glm::vec3 v;
    glm::vec3 normal;
    float D;
    glm::vec3 w;
    AABB bbox;
};





