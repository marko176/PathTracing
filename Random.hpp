#pragma once
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtx/intersect.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <random>
inline double random_double() {
    static std::uniform_real_distribution<double> distribution(0.0, 1.0);
    static thread_local std::mt19937 generator(std::random_device{}());
    return distribution(generator);
}

inline float random_float() {
    static std::uniform_real_distribution<float> distribution(0.0f, std::nextafterf(1.0f,0.0f));
    static thread_local std::mt19937 generator(std::random_device{}());
    return distribution(generator);
}

inline double random_double(double min, double max) {
    // Returns a random real in [min,max).
    //std::uniform_real_distribution(min,max);
    return min + (max-min)*random_double();
}

inline float random_float(float min, float max) {
    // Returns a random real in [min,max).
    //std::uniform_real_distribution(min,max);
    return min + (max-min)*random_float();
}

inline glm::vec3 random_unit_vector() {
    while (true) {
        glm::vec3 p = glm::vec3(random_double(-1,1),random_double(-1,1),random_double(-1,1));
        float lensq = p.x*p.x + p.y*p.y + p.z*p.z;
        if (lensq >= 1e-30f && lensq <= 1.0f){
            return glm::normalize(p);
        }
    }
}

inline glm::vec3 random_in_unit_disk() {
    while (true) {
        glm::vec3 p = glm::vec3(random_double(-1,1),random_double(-1,1),0);
        if (p.x*p.x + p.y*p.y <= 1.0f){
            return p;
        }
    }
}
inline glm::vec3 random_on_hemisphere(const glm::vec3& normal) {
    glm::vec3 vec = random_unit_vector();
    if(glm::dot(vec,normal)>0.0f){
        return vec;
    }else{
        return -vec;
    }
}

inline float schlick(float cosine, float ref_idx) {
    float r0 = (1-ref_idx) / (1+ref_idx);
    r0 = r0*r0;
    return r0 + (1-r0)*pow((1 - cosine),5.f);
}