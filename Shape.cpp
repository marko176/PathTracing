#include "Shape.hpp"
#include "Mesh.hpp"
bool SphereShape::Intersect(const Ray& ray, SurfaceInteraction& interaction, float max) const {
    glm::vec3 oc = ray.origin - center;
    float a = glm::dot(ray.dir, ray.dir);
    float b = glm::dot(oc, ray.dir);
    float c = glm::dot(oc, oc) - radius*radius;
    float discriminant = b*b - a*c;
    if (discriminant > 0) {
        float temp = (-b - std::sqrt(discriminant)) / a;
        if (temp < max && temp > shadowEpsilon) {
            interaction.t = temp;
            interaction.ns = glm::normalize(ray.at(temp) - center);
            interaction.n = interaction.ns;
            interaction.tangent = glm::vec3(0,0,0);
            interaction.p = ray.at(temp) + shadowEpsilon * interaction.n;
            interaction.uv = getSphereUV(interaction.n);
            return true;
        }
        temp = (-b + std::sqrt(discriminant)) / a;
        if (temp < max && temp > shadowEpsilon) {
            interaction.t = temp;
            interaction.ns = glm::normalize(ray.at(temp) - center);
            interaction.n = interaction.ns;
            interaction.tangent = glm::vec3(0,0,0);
            interaction.p = ray.at(temp) + shadowEpsilon * interaction.n;
            interaction.uv = getSphereUV(interaction.n);
            return true;
        }
    }
    return false;
}
bool SphereShape::IntersectPred(const Ray& ray, float max) const {
    glm::vec3 oc = ray.origin - center;
    float a = glm::dot(ray.dir, ray.dir);
    float b = glm::dot(oc, ray.dir);
    float c = glm::dot(oc, oc) - radius*radius;
    float discriminant = b*b - a*c;
    if (discriminant > 0) {
        float temp = (-b - std::sqrt(discriminant)) / a;
        if (temp < max && temp > shadowEpsilon) {
            return true;
        }
        temp = (-b + std::sqrt(discriminant)) / a;
        if (temp < max && temp > shadowEpsilon) {
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
    if(area * light_cosine == 0)return 0;
    return (dist_squared) / (light_cosine * area);
}

SurfaceInteraction SphereShape::Sample(const glm::vec2& u) const  {
    // sample a point uniformly on unit‐sphere
    float z    = 1.0f - 2.0f * u.x;
    float r    = std::sqrt(1.0f - z*z);
    float phi  = 2.0f * std::numbers::pi_v<float> * u.y;
    glm::vec3 dir{ r * std::cos(phi), r * std::sin(phi), z };
    glm::vec3 p = center + radius * dir;  
    return SurfaceInteraction{p, glm::normalize(p-center), getSphereUV(p)};      
}

inline bool IntersectRayTriangle(const glm::vec3& origin, const glm::vec3& dir,const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& v3, float& u,float& v, float& t) {
    glm::vec3 edge1 = v2 - v1;
    glm::vec3 edge2 = v3 - v1;
    glm::vec3 h = glm::cross(dir, edge2);
    float det = glm::dot(edge1,h);
    if(det > -std::numeric_limits<float>::epsilon() && det < std::numeric_limits<float>::epsilon())return false;
    float inv_det = 1.0f/det;
    glm::vec3 s = origin - v1;
    u = glm::dot(s,h) * inv_det;
    if(u < 0 || u > 1)return false;
    glm::vec3 q = glm::cross(s,edge1);
    v = glm::dot(dir,q) * inv_det;
    if(v < 0 || u + v > 1)return false;
    t = glm::dot(edge2, q) * inv_det;
    return t >= shadowEpsilon;
}


void sse_cross(__m128 result[3], const __m128 a[3], const __m128 b[3]){
    result[0] = _mm_fmsub_ps(a[1], b[2], _mm_mul_ps(b[1], a[2]));
    result[1] = _mm_fmsub_ps(a[2], b[0], _mm_mul_ps(b[2], a[0]));
    result[2] = _mm_fmsub_ps(a[0], b[1], _mm_mul_ps(b[0], a[1]));
}

__m128 sse_dot(const __m128 a[3], const __m128 b[3]){
    return _mm_fmadd_ps(a[2], b[2], _mm_fmadd_ps(a[1], b[1], _mm_mul_ps(a[0],b[0])));
}

void sse_sub(__m128 result[3], const __m128 a[3], const __m128 b[3]){
    result[0] = _mm_sub_ps(a[0],b[0]);
    result[1] = _mm_sub_ps(a[1],b[1]);
    result[2] = _mm_sub_ps(a[2],b[2]);
}


inline bool IntersectRay4Triangle(const glm::vec3& origin, const glm::vec3& dir,const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& v3, float& uu,float& vv, float& tt) {
    static const __m128 zeros = _mm_set1_ps(0);
    static const __m128 ones = _mm_set1_ps(1);
    static const __m128 minusOnes = _mm_set1_ps(-1);
    static const __m128 negativeEpsilon128 = _mm_set1_ps(-std::numeric_limits<float>::epsilon());
    static const __m128 positiveEpsilon128 = _mm_set1_ps(std::numeric_limits<float>::epsilon());
    __m128 mv1[3] = {_mm_set_ps(0,0,0,v1.x),_mm_set_ps(0,0,0,v1.y),_mm_set_ps(0,0,0,v1.z)};
    __m128 mv2[3] = {_mm_set_ps(0,0,0,v2.x),_mm_set_ps(0,0,0,v2.y),_mm_set_ps(0,0,0,v2.z)};
    __m128 mv3[3] = {_mm_set_ps(0,0,0,v3.x),_mm_set_ps(0,0,0,v3.y),_mm_set_ps(0,0,0,v3.z)};
    
    __m128 me1[3] = {_mm_sub_ps(mv2[0],mv1[0]),_mm_sub_ps(mv2[1],mv1[1]),_mm_sub_ps(mv2[2],mv1[2])};// glm::vec3 edge1 = v2 - v1;
    __m128 me2[3] = {_mm_sub_ps(mv3[0],mv1[0]),_mm_sub_ps(mv3[1],mv1[1]),_mm_sub_ps(mv3[2],mv1[2])};// glm::vec3 edge2 = v3 - v1;
    __m128 inactiveMask = _mm_set_ps(-1,-1,-1,0);

    __m128 mD4[3] = {_mm_set1_ps(dir.x),_mm_set1_ps(dir.y),_mm_set1_ps(dir.z)};
    __m128 mO4[3] = {_mm_set1_ps(origin.x),_mm_set1_ps(origin.y),_mm_set1_ps(origin.z)};
    __m128 mh[3];
    sse_cross(mh,mD4, me2);// glm::vec3 h = glm::cross(dir, edge2);
    
    __m128 mDet = sse_dot(me1,mh);// float det = glm::dot(edge1,h);
    __m128 minvDet = _mm_div_ps(ones,mDet);// float inv_det = 1.0f/det;
    __m128 ms[3];
    sse_sub(ms,mO4,mv1);// glm::vec3 s = origin - v1;

    __m128 u = _mm_mul_ps(minvDet,sse_dot(ms,mh));// u = glm::dot(s,h) * inv_det;
    __m128 mq[3];
    sse_cross(mq,ms,me1);// glm::vec3 q = glm::cross(s,edge1);
    __m128 v = _mm_mul_ps(minvDet,sse_dot(mD4,mq));// v = glm::dot(dir,q) * inv_det;

    __m128 t = _mm_mul_ps(minvDet,sse_dot(me2,mq));// t = glm::dot(edge2, q) * inv_det;

    //if(det > -std::numeric_limits<float>::epsilon() && det < std::numeric_limits<float>::epsilon())return false;
    __m128 failMask = _mm_and_ps(_mm_cmp_ps(mDet,negativeEpsilon128, _CMP_GT_OQ), _mm_cmp_ps(mDet,positiveEpsilon128, _CMP_LT_OQ));
    failMask = _mm_or_ps(failMask, _mm_cmp_ps(u,zeros,_CMP_LT_OQ));//u < 0
    failMask = _mm_or_ps(failMask, _mm_cmp_ps(v,zeros,_CMP_LT_OQ));//v < 0
    failMask = _mm_or_ps(failMask, _mm_cmp_ps(_mm_add_ps(u,v),ones,_CMP_GT_OQ));//u + v > 1
    //t > max here
    failMask = _mm_or_ps(failMask, _mm_cmp_ps(t,zeros,_CMP_LT_OQ));
    failMask = _mm_or_ps(failMask,inactiveMask);

    __m128 tfinal = _mm_blendv_ps(t, minusOnes, failMask);

    int mask = _mm_movemask_ps(tfinal);
    if(mask != 0b1111){
        bool intersected = false;
        float* ptr = (float*)&tfinal;
        float* uptr = (float*)&u;
        float* vptr = (float*)&v;
        for(int i = 0;i<4;i++){
            if(ptr[i]>0 && ptr[i] < tt){
                tt = ptr[i];
                uu = uptr[i];
                vv = vptr[i];
                intersected = true;
            }
            
        }
        return intersected;
    }

    return false;
}



bool TriangleShape::Intersect(const Ray& ray, SurfaceInteraction& interaction,float max) const {
    
    glm::vec2 baryPos;
    float t = std::numeric_limits<float>::infinity();
    const Mesh* mesh = meshList[MeshIndex];
    int index0 = mesh->indices[TriIndex*3+0];
    int index1 = mesh->indices[TriIndex*3+1];
    int index2 = mesh->indices[TriIndex*3+2];

    bool hit_triangle = glm::intersectRayTriangle(ray.origin,ray.dir,mesh->vertices[index0]
                                                                    ,mesh->vertices[index1]
                                                                    ,mesh->vertices[index2],baryPos,t);
    //bool hit_triangle = IntersectRayTriangle(ray.origin,ray.dir,mesh->vertices[index0],mesh->vertices[index1],mesh->vertices[index2],baryPos.x,baryPos.y,t);
    //bool hit_triangle = IntersectRay4Triangle(ray.origin,ray.dir,mesh->vertices[index0],mesh->vertices[index1],mesh->vertices[index2],baryPos.x,baryPos.y,t);
    //if we dont use glm t<shadowEPsilon not needed
    if(!hit_triangle || t > max || t < shadowEpsilon)return false;
        
    float u = baryPos.x;
    float v = baryPos.y;
    float w = 1.0f - u - v;
    
    glm::vec2 uv =  u * mesh->texCoords[index1] +
                    v * mesh->texCoords[index2]+
                    w * mesh->texCoords[index0];

    if(random_float()>mesh->material->Alpha(uv.x,uv.y))return 0;
    glm::vec3 norm_normal = glm::normalize( u * mesh->normals[index1] +
                                            v * mesh->normals[index2] +
                                            w * mesh->normals[index0]);
                                            
    glm::vec3 e1 = mesh->vertices[index1] - mesh->vertices[index0];
    glm::vec3 e2 = mesh->vertices[index2] - mesh->vertices[index0];
    glm::vec3 N = glm::normalize(glm::cross(e1,e2));
    interaction.n = N;
    if(glm::dot(N,norm_normal)<0){
        norm_normal = -norm_normal;
    }
    if(glm::dot(ray.dir,N)>0.0f){
        //norm_normal = -norm_normal;//WRONG BUT TIR DOESNT WORK
        N*=-1;
    }
    interaction.t = t;
    interaction.ns = norm_normal;  
    interaction.uv = uv;
    interaction.p = ray.at(t) + shadowEpsilon * N;
    interaction.AreaLight = nullptr;
    interaction.mat = mesh->material;
    if(!mesh->tangents.empty()){
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

    return true;
}
bool TriangleShape::IntersectPred(const Ray& ray, float max) const {
    
    glm::vec2 baryPos;
    float t = std::numeric_limits<float>::infinity();
    const Mesh* mesh = meshList[MeshIndex];
    int index0 = mesh->indices[TriIndex*3+0];
    int index1 = mesh->indices[TriIndex*3+1];
    int index2 = mesh->indices[TriIndex*3+2];
    bool hit_triangle = glm::intersectRayTriangle(ray.origin,ray.dir,mesh->vertices[index0]
                                                                    ,mesh->vertices[index1]
                                                                    ,mesh->vertices[index2],baryPos,t);
    if(!hit_triangle || t > max || t < shadowEpsilon)return false;
        
    float u = baryPos.x;
    float v = baryPos.y;
    float w = 1.0f - u - v;
    
  
    glm::vec2 uv =  u * mesh->texCoords[index1] +
                    v * mesh->texCoords[index2]+
                    w * mesh->texCoords[index0];

   
    return random_float()<=mesh->material->Alpha(uv.x,uv.y);
}

AABB TriangleShape::BoundingBox() const {
    const Mesh* mesh = meshList[MeshIndex];
    AABB bbox(mesh->vertices[mesh->indices[TriIndex*3+0]]);
    bbox.Expand(mesh->vertices[mesh->indices[TriIndex*3+1]]);
    bbox.Expand(mesh->vertices[mesh->indices[TriIndex*3+2]]);
    return bbox;
}
SurfaceInteraction TriangleShape::Sample(const glm::vec2& u) const {
    float w = 1.0f - u.x - u.y;
    Mesh* mesh = meshList[MeshIndex];

    int index0 = mesh->indices[TriIndex*3+0];
    int index1 = mesh->indices[TriIndex*3+1];
    int index2 = mesh->indices[TriIndex*3+2];

    glm::vec3 e1 = mesh->vertices[index1] - mesh->vertices[index0];
    glm::vec3 e2 = mesh->vertices[index2] - mesh->vertices[index0];
    glm::vec3 n = glm::normalize(glm::cross(e1,e2));
    if(n.x != n.x){
        n = {0,0,0};
    }
    


    //glm::vec3 p = u.x * mesh->normals[mesh->indices[index1]] + u.y * mesh->normals[mesh->indices[index2]] + w * mesh->normals[mesh->indices[index0]];
    
    glm::vec3 p = u.x * mesh->vertices[index1] + u.y * mesh->vertices[index2] + w * mesh->vertices[index0];
    glm::vec2 uv =  u.x * mesh->texCoords[index1] +
                    u.y * mesh->texCoords[index2]+
                    w * mesh->texCoords[index0];
    return SurfaceInteraction{p,n,uv};
}
float TriangleShape::Area() const {   
    Mesh* mesh = meshList[MeshIndex];                                                   
    return glm::length(glm::cross(mesh->vertices[mesh->indices[TriIndex*3+0]] - mesh->vertices[mesh->indices[TriIndex*3+2]],mesh->vertices[mesh->indices[TriIndex*3+1]] - mesh->vertices[mesh->indices[TriIndex*3+2]]))*0.5f;
}

float TriangleShape::PDF(const GeometricInteraction& interaction) const {
    float area = Area();
    if(area == 0 || interaction.n.x != interaction.n.x)return 0;
    return 1.0f/area; //need to sample only visible part! -> put into another function for eg sampleCone, PDFCone ?
}
float TriangleShape::PDF(const GeometricInteraction& interaction,const Ray& ray) const {
    
    glm::vec3 to_shape = interaction.p - ray.origin;
    float dist_squared = glm::dot(to_shape,to_shape);
    float light_cosine = std::abs(glm::dot(-ray.dir,interaction.n));
    float area = Area();
    if(area == 0 || light_cosine == 0 || interaction.n.x != interaction.n.x)return 0;
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
    if(t < shadowEpsilon || t > max)return false;
    glm::vec3 planar_hit = ray.at(t) - Q;
    float alpha = glm::dot(w,glm::cross(planar_hit,v));
    float beta = glm::dot(w,glm::cross(u,planar_hit));
    if(!is_interior(alpha,beta))return false;
    interaction.uv = {alpha,beta};
    interaction.t = t;
    interaction.ns = norm_normal;
    interaction.n = normal;
    interaction.tangent = glm::vec3(0,0,0);
    interaction.p = ray.at(t) + shadowEpsilon * norm_normal;//was norm_normal
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
    if(t < shadowEpsilon || t > max)return false;
    glm::vec3 planar_hit = ray.at(t) - Q;
    float alpha = glm::dot(w,glm::cross(planar_hit,v));
    float beta = glm::dot(w,glm::cross(u,planar_hit));
    return is_interior(alpha,beta);
}