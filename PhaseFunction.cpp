#include "PhaseFunction.hpp"
#include "Onb.hpp"
float HenyeyGreenstein::PDF(const glm::vec3& in, const glm::vec3& out) const {
    return phaseHG<float>(glm::dot(in,out),g);
}


float HenyeyGreenstein::Sample(const glm::vec3& in, glm::vec3& out,const glm::vec2& u) const {
    float cosTheta;
    if (std::abs(g) < 1e-3)
        cosTheta = 1 - 2 * u[0];
    else {
        float sqrTerm = (1 - g * g) / (1 - g + 2 * g * u[0]);
        cosTheta = (1 + g * g - sqrTerm * sqrTerm) / (2 * g);
    }

    float sinTheta = std::sqrt(std::max<float>(0,
                                           1 - cosTheta * cosTheta));
    float phi = 2 * std::numbers::pi_v<float> * u[1];
    float x = sinTheta * std::cos(phi);
    float y = sinTheta * std::sin(phi);
    float z = cosTheta; // local +z is "forward" (incoming)

    /*
    glm::vec3 forward = in; // incoming direction at interaction
    glm::vec3 up = (std::abs(forward.z) < 0.9999f) ? glm::vec3(1.f,0.f,0.f) : glm::vec3(0.f,1.f,0.f);
    glm::vec3 right = glm::normalize(glm::cross(up, forward));
    glm::vec3 newUp = glm::cross(forward, right);
    */

    onb basis(in);


    out = glm::normalize(basis.toWorld({x,y,z}));

    return phaseHG<float>(cosTheta, g);
}