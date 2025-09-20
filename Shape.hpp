#pragma once
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtx/intersect.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <memory>
#include "Ray.hpp"
#include "Hit_record.hpp"
#include "Material.hpp"
#include <mutex>

class Mesh;
class Shape {
public:
    virtual AABB Bounding_box() const = 0;
    virtual bool IntersectPred(const Ray& ray, float max = 1e30f) const = 0;
    virtual bool Intersect(const Ray& ray, SurfaceInteraction& interaction, float max = 1e30f) const = 0;
    virtual float Area() const = 0;
    virtual GeometricInteraction Sample(const glm::vec2& u) const = 0;
    virtual float PDF(const GeometricInteraction& interaction) const = 0;
    virtual float PDF(const GeometricInteraction& interaction,const Ray& ray) const = 0; //should just give normal info (shape sample context, maybe get that from sample?) 
};

class SphereShape : public Shape{ 
public:
    SphereShape(const glm::vec3& center, float r) : center(center), bbox{}, radius(r) {
        glm::vec3 rvec = glm::vec3(radius);
        bbox.expand(center-rvec);
        bbox.expand(center+rvec);
    }
    bool Intersect(const Ray& ray, SurfaceInteraction& interaction, float max = 1e30f) const override ;

    bool IntersectPred(const Ray& ray, float max = 1e30f) const override ;

    AABB Bounding_box() const override{
        return bbox;
    }

    static glm::vec2 getSphereUV(const glm::vec3& p){
        float theta = std::acos(std::clamp(p.y,-1.0f,1.0f));
        float phi = std::atan2f(p.z,p.x);
        if(phi < 0) phi += 2.0f * std::numbers::pi_v<float>;
        float u = std::numbers::inv_pi_v<float> * phi * 0.5f;
        float v = std::numbers::inv_pi_v<float> * theta;
        return {u,v};
    }


    GeometricInteraction Sample(const glm::vec2& u) const override ;

    float Area() const override ;

    float PDF(const GeometricInteraction& interaction) const override ;

    float PDF(const GeometricInteraction& interaction,const Ray& ray) const override ;

private:
    glm::vec3 center;
    AABB bbox;
    float radius;
};

class TriangleShape : public Shape {
public:
    TriangleShape(uint32_t meshIndex, uint32_t TriIndex) : MeshIndex(meshIndex), TriIndex(TriIndex) {}

    bool Intersect(const Ray& ray, SurfaceInteraction& interaction, float max = 1e30f) const override;

    bool IntersectPred(const Ray& ray, float max = 1e30f) const override;

    AABB Bounding_box() const override;

    GeometricInteraction Sample(const glm::vec2& u) const override;

    float Area() const override;

    float PDF(const GeometricInteraction& interaction) const override;

    float PDF(const GeometricInteraction& interaction,const Ray& ray) const override;

    //maybe move this to resource manager?
    [[nodiscard]] static std::size_t addMesh(Mesh* mesh) {
        const std::lock_guard<std::mutex> ml(meshListLock);
        for(int i = 0;i<meshList.size();i++){
            if(meshList[i]==nullptr || meshList[i] == mesh){
                meshList[i]=mesh;
                return i;
            }
        }
        meshList.push_back(mesh);
        return meshList.size()-1;
    }
    
    static void removeMesh(Mesh* mesh) {
        const std::lock_guard<std::mutex> ml(meshListLock);
        for(int i = 0;i<meshList.size();i++){
            if(meshList[i]==mesh){
                meshList[i]=nullptr;
            }
        }
    }
private:
    uint32_t MeshIndex;
    uint32_t TriIndex;
    static inline std::vector<Mesh*> meshList;
    static inline std::mutex meshListLock;
};

class QuadShape : public Shape {
public:
    QuadShape(const glm::vec3& Q, const glm::vec3& u, const glm::vec3& v) : Q(Q), u(u), v(v) , bbox{} {
        glm::vec3 n = glm::cross(u,v);
        normal = glm::normalize(n);
        D = glm::dot(normal,Q);
        bbox.expand(Q);
        bbox.expand(Q + u + v);
        bbox.expand(Q + u);
        bbox.expand(Q + v);
        w = n / glm::dot(n,n);
    }

    bool Intersect(const Ray& ray, SurfaceInteraction& interaction, float max = 1e30f) const override ;

    bool IntersectPred(const Ray& ray, float max = 1e30f) const override ;

    AABB Bounding_box() const override{
        return bbox;
    }

    GeometricInteraction Sample(const glm::vec2& uv) const override{
        return GeometricInteraction{Q + uv.x * u + uv.y * v, normal, uv};
    }

    float Area() const override{                                                      
        return glm::length(glm::cross(u,v));
    }

    float PDF(const GeometricInteraction& interaction) const override {
        return 1.0f/Area();
    }

    float PDF(const GeometricInteraction& interaction,const Ray& ray) const override {
        glm::vec3 to_shape = interaction.p - ray.origin;
        float dist_squared = glm::dot(to_shape,to_shape);
        float light_cosine = std::abs(glm::dot(-ray.dir,interaction.n));
        float area = Area();
        if(area == 0)return 0;
        return dist_squared / (light_cosine * area);
    }
private:

    static bool is_interior(float a, float b, GeometricInteraction& interaction) {
        if (a<0 || a>1 || b<0 || b>1)
            return false;
        interaction.uv = {a,b};
        return true;
    }

    static bool is_interior(float a, float b) {
        return a>=0 && a<=1 && b>=0 && b<=1;
    }
    glm::vec3 Q;
    glm::vec3 u,v;
    glm::vec3 normal;
    float D;
    glm::vec3 w;
    AABB bbox;
};





