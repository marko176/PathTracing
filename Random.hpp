#pragma once
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtx/intersect.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <random>
#include <numbers>
inline double random_double() {
    static thread_local std::uniform_real_distribution<double> distribution(0.0, 1.0);
    static thread_local std::mt19937 generator(std::random_device{}());
    return distribution(generator);
}

inline float random_float() {
    static thread_local std::uniform_real_distribution<float> distribution(0.0f, std::nextafterf(1.0f,0.0f));
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


inline glm::vec2 inUnitDisk(const glm::vec2& uv) {
    float r = std::sqrt(uv.x);
    float theta = 2 * std::numbers::pi_v<float> * uv.y;
    return r*glm::vec2{std::cos(theta),std::sin(theta)};
}


