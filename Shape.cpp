#include "Shape.hpp"
#include "Mesh.hpp"
bool SphereShape::Intersect(const Ray& ray, SurfaceInteraction& interaction, float max) const {
    glm::vec3 oc = ray.origin - center;
    float a = dot(ray.dir, ray.dir);
    float b = dot(oc, ray.dir);
    float c = dot(oc, oc) - radius*radius;
    float discriminant = b*b - a*c;
    if (discriminant > 0) {
        float temp = (-b - std::sqrt(discriminant)) / a;
        if (temp < max && temp > 0.0001f) {
            interaction.t = temp;
            interaction.ns = (ray.at(temp) - center) / radius;
            interaction.n = interaction.ns;
            interaction.tangent = glm::vec3(0,0,0);
            interaction.p = ray.at(temp) + 0.0001f * interaction.n;
            interaction.uv = getSphereUV(interaction.n);
            return true;
        }
        temp = (-b + std::sqrt(discriminant)) / a;
        if (temp < max && temp > 0.0001f) {
            interaction.t = temp;
            interaction.ns = (ray.at(temp) - center) / radius;
            interaction.n = interaction.ns;
            interaction.tangent = glm::vec3(0,0,0);
            interaction.p = ray.at(temp) + 0.0001f * interaction.n;
            interaction.uv = getSphereUV(interaction.n);
            return true;
        }
    }
    return false;
}
bool SphereShape::IntersectPred(const Ray& ray, float max) const {
    glm::vec3 oc = ray.origin - center;
    float a = dot(ray.dir, ray.dir);
    float b = dot(oc, ray.dir);
    float c = dot(oc, oc) - radius*radius;
    float discriminant = b*b - a*c;
    if (discriminant > 0) {
        float temp = (-b - std::sqrt(discriminant)) / a;
        if (temp < max && temp > 0.0001f) {
            return true;
        }
        temp = (-b + std::sqrt(discriminant)) / a;
        if (temp < max && temp > 0.0001f) {
            return true;
        }
    }
    return false;
}

float SphereShape::Area() const{
    return 4.0f*std::numbers::pi_v<float>*radius*radius;
}
float SphereShape::PDF(const GeometricInteraction& interaction) const {
    return 1.0f/Area();
}
float SphereShape::PDF(const GeometricInteraction& interaction,const Ray& ray) const {
    glm::vec3 to_shape = interaction.p - ray.origin;
    float dist_squared = glm::dot(to_shape,to_shape);
    float light_cosine = std::abs(glm::dot(-ray.dir,interaction.n));
    float area = Area();
    if(area == 0)return 0;
    return (dist_squared) / (light_cosine * area);
}

GeometricInteraction SphereShape::Sample(const glm::vec2& u) const  {
    // sample a point uniformly on unit‐sphere
    float z    = 1.0f - 2.0f * u.x;
    float r    = std::sqrt(1.0f - z*z);
    float phi  = 2.0f * std::numbers::pi_v<float> * u.y;
    glm::vec3 dir{ r * std::cos(phi), r * std::sin(phi), z };
    glm::vec3 p = center + radius * dir;  
    return GeometricInteraction{p, (p-center)/ radius};      
}



bool TriangleShape::Intersect(const Ray& ray, SurfaceInteraction& interaction,float max) const {
    
    glm::vec2 baryPos;
    float t = -1;
    Mesh* mesh = meshList[MeshIndex];
    int index0 = mesh->indices[TriIndex*3+0];
    int index1 = mesh->indices[TriIndex*3+1];
    int index2 = mesh->indices[TriIndex*3+2];
    bool hit_triangle = glm::intersectRayTriangle(ray.origin,ray.dir,mesh->vertices[index0]
                                                                    ,mesh->vertices[index1]
                                                                    ,mesh->vertices[index2],baryPos,t);
    if(!hit_triangle || t > max || t < 0.0001f)return false;
        
    float u = baryPos.x;
    float v = baryPos.y;
    float w = 1.0f - u - v;
    
    glm::vec2 uv =  u * mesh->texCoords[index1] +
                    v * mesh->texCoords[index2]+
                    w * mesh->texCoords[index0];

    if(random_float()>mesh->material->Alpha(u,v))return 0;
    glm::vec3 norm_normal = glm::normalize( u * mesh->normals[index1] +
                                            v * mesh->normals[index2] +
                                            w * mesh->normals[index0]);
                                            
    glm::vec3 e1 = mesh->vertices[index1] - mesh->vertices[index0];
    glm::vec3 e2 = mesh->vertices[index2] - mesh->vertices[index0];
    glm::vec3 N = glm::normalize(glm::cross(e1,e2));
    interaction.n = N;
    if(glm::dot(ray.dir,N)>0.0f){
        norm_normal = -norm_normal;
        N*=-1;
    }
    interaction.t = t;
    interaction.ns = norm_normal;  
    interaction.uv = uv;

    //if angle small do just ray.at(t- eps) ? 
    interaction.p = ray.at(t) + 0.0001f * N;//was 0.0001 * N should be 0.0001f for dragon
    //interaction.n = N;
    interaction.AreaLight = nullptr;//(hittable*)&primitives[TriIndex];
    interaction.mat = mesh->material;
    if(!mesh->tangents.empty() /*&& !mesh->bitangents.empty() */){
        glm::vec3 tangent =  u * mesh->tangents[index1] +
                            v * mesh->tangents[index2]+
                            w * mesh->tangents[index0];
        interaction.tangent = glm::normalize(tangent);
        glm::vec3 bitangent =  u * mesh->bitangents[index1] +
                                v * mesh->bitangents[index2]+
                                w * mesh->bitangents[index0];
        interaction.bitangent = glm::normalize(bitangent);
        //we should jsut get pointer to texture then not store tangent and bitangent in interaction!!! pbrt does that -> same for albedo?
        interaction.ns = interaction.mat->sample_normalMap(interaction);
    }else{
        interaction.tangent = glm::vec3(0,0,0);
    }
    //if(glm::dot(rec.normal,ray.dir)>=0)rec.normal = norm_normal;
    return true;
}
bool TriangleShape::IntersectPred(const Ray& ray, float max) const {
    
    glm::vec2 baryPos;
    float t = -1;
    Mesh* mesh = meshList[MeshIndex];
    int index0 = mesh->indices[TriIndex*3+0];
    int index1 = mesh->indices[TriIndex*3+1];
    int index2 = mesh->indices[TriIndex*3+2];
    bool hit_triangle = glm::intersectRayTriangle(ray.origin,ray.dir,mesh->vertices[index0]
                                                                    ,mesh->vertices[index1]
                                                                    ,mesh->vertices[index2],baryPos,t);
    if(!hit_triangle || t > max || t < 0.0001f)return false;
        
    float u = baryPos.x;
    float v = baryPos.y;
    float w = 1.0f - u - v;
    
    //mesh->normal.at(u,v);
    glm::vec3 norm_normal = glm::normalize( u * mesh->normals[index1] +
                                            v * mesh->normals[index2] +
                                            w * mesh->normals[index0]);
                                            
    //if(mesh->normals[mesh->indices[TriIndex*3]] != mesh->normals[index1])std::cout<<"FAIL";
    glm::vec3 e1 = mesh->vertices[index1] - mesh->vertices[index0];
    glm::vec3 e2 = mesh->vertices[index2] - mesh->vertices[index0];
    glm::vec3 N = glm::normalize(glm::cross(e1,e2));
    if(glm::dot(ray.dir,N)>0.0f){
        norm_normal = -norm_normal;
        N*=-1;
    }
  
    glm::vec2 uv =  u * mesh->texCoords[index1] +
                    v * mesh->texCoords[index2]+
                    w * mesh->texCoords[index0];

    //if(random_float()>mesh->material->Alpha(u,v))return false;
    return random_float()<=mesh->material->Alpha(u,v);
}

AABB TriangleShape::BoundingBox() const {
    Mesh* mesh = meshList[MeshIndex];
    AABB bbox(mesh->vertices[mesh->indices[TriIndex*3+0]]);
    bbox.Expand(mesh->vertices[mesh->indices[TriIndex*3+1]]);
    bbox.Expand(mesh->vertices[mesh->indices[TriIndex*3+2]]);
    return bbox;
}
GeometricInteraction TriangleShape::Sample(const glm::vec2& u) const {
    float w = 1.0f - u.x - u.y;
    Mesh* mesh = meshList[MeshIndex];

    glm::vec3 e1 = mesh->vertices[TriIndex*3+1] - mesh->vertices[TriIndex*3+0];
    glm::vec3 e2 = mesh->vertices[TriIndex*3+2] - mesh->vertices[TriIndex*3+0];
    glm::vec3 n = glm::normalize(glm::cross(e1,e2));
    glm::vec3 p = u.x * mesh->normals[mesh->indices[TriIndex*3+1]] + u.y * mesh->normals[mesh->indices[TriIndex*3+2]] + w * mesh->normals[mesh->indices[TriIndex*3+0]];
    return GeometricInteraction{p,n};
}
float TriangleShape::Area() const {   
    Mesh* mesh = meshList[MeshIndex];                                                   
    return glm::length(glm::cross(mesh->vertices[mesh->indices[TriIndex*3+0]] - mesh->vertices[mesh->indices[TriIndex*3+2]],mesh->vertices[mesh->indices[TriIndex*3+1]] - mesh->vertices[mesh->indices[TriIndex*3+2]]))*0.5f;
}

float TriangleShape::PDF(const GeometricInteraction& interaction) const {
    return 1.0f/Area(); //need to sample only visible part! -> put into another function for eg sampleCone, PDFCone ?
}
float TriangleShape::PDF(const GeometricInteraction& interaction,const Ray& ray) const {
    glm::vec3 to_shape = interaction.p - ray.origin;
    float dist_squared = glm::dot(to_shape,to_shape);
    float light_cosine = std::abs(glm::dot(-ray.dir,interaction.n));
    float area = Area();
    if(area == 0 || light_cosine == 0)return 0;
    return (dist_squared) / (light_cosine * area);
}


bool QuadShape::Intersect(const Ray& ray, SurfaceInteraction& interaction,float max) const {
    glm::vec3 norm_normal = normal;
    float DD = D;
    if(glm::dot(ray.dir,normal)>0){
        norm_normal = -normal;
        DD = -D;
    }
    float denom = glm::dot(norm_normal,ray.dir);
    if(std::fabs(denom) < 1e-8f)return false;
    float t = (DD - glm::dot(norm_normal,ray.origin)) / denom;
    if(t < 0.0001f || t > max)return false;
    glm::vec3 planar_hit = ray.at(t) - Q;
    float alpha = glm::dot(w,glm::cross(planar_hit,v));
    float beta = glm::dot(w,glm::cross(u,planar_hit));
    if(!is_interior(alpha,beta))return false;
    interaction.uv = {alpha,beta};
    interaction.t = t;
    interaction.ns = norm_normal;
    interaction.n = normal;
    interaction.tangent = glm::vec3(0,0,0);
    interaction.p = ray.at(t) + 0.0001f * norm_normal;//was norm_normal
    return true;
}
bool QuadShape::IntersectPred(const Ray& ray, float max) const {
    glm::vec3 norm_normal = normal;
    float DD = D;
    if(glm::dot(ray.dir,normal)>0){
        norm_normal = -normal;
        DD = -D;
    }
    float denom = glm::dot(norm_normal,ray.dir);
    if(std::fabs(denom) < 1e-8f)return false;
    float t = (DD - glm::dot(norm_normal,ray.origin)) / denom;
    if(t < 0.0001f || t > max)return false;
    glm::vec3 planar_hit = ray.at(t) - Q;
    float alpha = glm::dot(w,glm::cross(planar_hit,v));
    float beta = glm::dot(w,glm::cross(u,planar_hit));
    return is_interior(alpha,beta);
}