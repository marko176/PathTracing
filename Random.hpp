#pragma once
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <random>
#include <numbers>
#include "pcg-cpp/pcg_random.hpp"



inline double random_double(){
    static thread_local std::uniform_real_distribution<double> distribution(0.0, 1.0);
    static thread_local pcg64 generator(std::random_device {}());
    return distribution(generator);
}

inline float random_float(){
    static thread_local std::uniform_real_distribution<float> distribution(0.0f, 1.0f);
    static thread_local pcg32 generator(std::random_device {}());
    return distribution(generator);
}

inline std::array<float, 8> RandomFloatx8(){
    alignas(32) std::array<float, 8> tmp;
#if defined(__SSE__) || defined(_M_AMD64) || defined(_M_X64)
    static thread_local pcg32 generator(std::random_device {}());
    static const __m256 scale = _mm256_set1_ps(1.0f / 0xFFFFFFFFu);

    const __m256 ints = _mm256_set_ps(static_cast<float>(generator()),
        static_cast<float>(generator()),
        static_cast<float>(generator()),
        static_cast<float>(generator()),
        static_cast<float>(generator()),
        static_cast<float>(generator()),
        static_cast<float>(generator()),
        static_cast<float>(generator()));
    const __m256 floats = _mm256_mul_ps(ints, scale);
    _mm256_store_ps(tmp.data(), floats);
#else
    for(int i = 0;i<8;++i){
        tmp[i]=random_float();
    }
#endif
    return tmp;
}


inline double random_double(double min, double max){
    // Returns a random real in [min,max).
    //std::uniform_real_distribution(min,max);
    return min + (max - min) * random_double();
}

inline float random_float(float min, float max){
    // Returns a random real in [min,max).
    //std::uniform_real_distribution(min,max);
    return min + (max - min) * random_float();
}


inline glm::vec2 inUnitDisk(const glm::vec2& uv){
    float r = std::sqrt(uv.x);
    float theta = 2 * std::numbers::pi_v<float> *uv.y;
    return r * glm::vec2 { std::cos(theta),std::sin(theta) };
}


