#pragma once
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtx/intersect.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <memory>
#include "Ray.hpp"
#include "Hit_record.hpp"
#include "Material.hpp"

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
    SphereShape(const glm::vec3& center,const std::shared_ptr<Material>& m, float r) : center(center), bbox{}, radius(r) {
        glm::vec3 rvec = glm::vec3(radius);
        bbox.expand(center-rvec);
        bbox.expand(center+rvec);
    }
    virtual bool Intersect(const Ray& ray, SurfaceInteraction& interaction, float max = 1e30f) const override ;

    virtual bool IntersectPred(const Ray& ray, float max = 1e30f) const override ;

    virtual AABB Bounding_box() const override{
        return bbox;
    }

    static void get_sphere_uv(const glm::vec3& p, float& u,float& v){
        float theta = std::acos(-p.y);
        float phi = std::atan2f(-p.z,p.x) + std::numbers::pi_v<float>;
        u = phi / (2*std::numbers::pi_v<float>);
        v = theta / std::numbers::pi_v<float>;
    }


    virtual GeometricInteraction Sample(const glm::vec2& u) const override ;

    virtual float Area() const override ;

    virtual float PDF(const GeometricInteraction& interaction) const override ;

    virtual float PDF(const GeometricInteraction& interaction,const Ray& ray) const override ;


private:
    glm::vec3 center;
    AABB bbox;
    float radius;
};

class TriangleShape : public Shape {

public:
    TriangleShape(uint32_t meshIndex, uint32_t TriIndex) : MeshIndex(meshIndex), TriIndex(TriIndex) {

    }

    virtual bool Intersect(const Ray& ray, SurfaceInteraction& interaction, float max = 1e30f) const override ;

    virtual bool IntersectPred(const Ray& ray, float max = 1e30f) const override ;

    virtual AABB Bounding_box() const override;

    virtual GeometricInteraction Sample(const glm::vec2& u) const override;

    virtual float Area() const override;

    virtual float PDF(const GeometricInteraction& interaction) const override ;

    virtual float PDF(const GeometricInteraction& interaction,const Ray& ray) const override ;

    static void addMesh(Mesh* mesh) {
        for(int i = 0;i<meshList.size();i++){
            if(meshList[i]==nullptr){
                meshList[i]=mesh;
                break;
            }
        }
        meshList.push_back(mesh);
    }
    
    static void removeMesh(Mesh* mesh) {
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

    virtual bool Intersect(const Ray& ray, SurfaceInteraction& interaction, float max = 1e30f) const override ;

    virtual bool IntersectPred(const Ray& ray, float max = 1e30f) const override ;

    virtual AABB Bounding_box() const override{
        return bbox;
    }

    virtual GeometricInteraction Sample(const glm::vec2& uv) const override{
        return GeometricInteraction{Q + uv.x * u + uv.y * v, normal, uv};
    }

    virtual float Area() const override{                                                      
        return glm::length(glm::cross(u,v));;
    }

    virtual float PDF(const GeometricInteraction& interaction) const override {
        return 1.0f/Area(); //need to sample only visible part! -> put into another function for eg sampleCone, PDFCone ?
    }

    virtual float PDF(const GeometricInteraction& interaction,const Ray& ray) const override {
        glm::vec3 to_shape = interaction.p - ray.origin;
        float dist_squared = glm::dot(to_shape,to_shape);
        float light_cosine = std::abs(glm::dot(-ray.dir,interaction.n));
        float area = Area();
        if(area == 0)return 0;
        return (dist_squared) / (light_cosine * area);
    }
private:

    static bool is_interior(float a, float b, GeometricInteraction& interaction) {
    
        // Given the hit point in plane coordinates, return false if it is outside the
        // primitive, otherwise set the hit record UV coordinates and return true.

        if (a<0 || a>1 || b<0 || b>1)
            return false;

        interaction.uv = {a,b};
        return true;
    }

    static bool is_interior(float a, float b) {
    
        // Given the hit point in plane coordinates, return false if it is outside the
        // primitive, otherwise set the hit record UV coordinates and return true.

        if (a<0 || a>1 || b<0 || b>1)
            return false;


        return true;
    }
    glm::vec3 Q;
    glm::vec3 u,v;
    glm::vec3 normal;
    float D;
    glm::vec3 w;
    AABB bbox;

};





