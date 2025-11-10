#include "PhaseFunction.hpp"
#include "Onb.hpp"

float HenyeyGreenstein::PDF(const glm::vec3& in, const glm::vec3& out) const{
    return phaseHG<float>(glm::dot(in, out), g);
}


float HenyeyGreenstein::Sample(const glm::vec3& in, glm::vec3& out, const glm::vec2& u) const{
    float cosTheta;
    if(std::abs(g) < 1e-3)
        cosTheta = 1 - 2 * u[0];
    else{
        float sqr = (1 - g * g) / (1 - g + 2 * g * u[0]);
        cosTheta = (1 + g * g - sqr * sqr) / (2 * g);
    }

    float sinTheta = std::sqrt(std::max<float>(0, 1 - cosTheta * cosTheta));
    float phi = 2 * std::numbers::pi_v<float> *u[1];
    float x = sinTheta * std::cos(phi);
    float y = sinTheta * std::sin(phi);
    float z = cosTheta; // local +z is "forward" (incoming)

    out = glm::normalize(onb(in).toWorld({ x,y,z }));

    return phaseHG<float>(cosTheta, g);
}