#pragma once
#include "Interaction.hpp"
#include <numbers>
template <std::floating_point T>
inline T phaseHG(T cosTheta, T g) {
    T denom = 1 + g * g + 2 * g * cosTheta;
    return static_cast<T>(0.25) * std::numbers::inv_pi_v<T> * (static_cast<T>(1) - g * g) / (denom * std::sqrt(denom));
}

class PhaseFunction {
public:
    virtual ~PhaseFunction() = default;
    virtual float PDF(const glm::vec3& in, const glm::vec3& out) const = 0;
    virtual float Sample(const glm::vec3& in, glm::vec3& out, const glm::vec2& u) const = 0;
};

class HenyeyGreenstein : public PhaseFunction {
public:
    virtual ~HenyeyGreenstein() = default;
    HenyeyGreenstein(float g) : g(g) {}

    float PDF(const glm::vec3& in, const glm::vec3& out) const override;
    float Sample(const glm::vec3& in, glm::vec3& out, const glm::vec2& u) const override;
private:
    float g;
};