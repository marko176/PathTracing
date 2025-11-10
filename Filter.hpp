#pragma once
#include "Interaction.hpp"
#include <concepts>
#include <algorithm>
#include <numbers>
template <std::floating_point T>
inline constexpr T Gaussian(T x, T sigma){
    //return 1.0 / std::sqrt( 2 * std::numbers::pi_v<T> * sigma * sigma) * std::exp(-(x*x) / (2 * sigma * sigma));
    return std::numbers::inv_sqrtpi_v<T> / (sigma * std::numbers::sqrt2_v<T>) * std::exp(-(x * x) / (2 * sigma * sigma));
}

template <std::floating_point T>
inline constexpr T GaussianIntegral(T x0, T x1, T sigma){
    T sigmaRoot2 = sigma * std::numbers::sqrt2_v<T>;
    return static_cast<T>(0.5) * (std::erf(-x0 / sigmaRoot2) - std::erf(-x1 / sigmaRoot2));
}

template <std::floating_point T>
inline constexpr T Sinc(T x){
    if(static_cast<T>(1) - x * x == static_cast<T>(1))return static_cast<T>(1);
    return std::sin(std::numbers::pi_v<T> *x) / (std::numbers::pi_v<T> *x);
}

template <std::floating_point T>
inline constexpr T WindowedSinc(T x, T radius, T tau){
    if(std::abs(x) > radius)return 0;
    return Sinc(x) * Sinc(x / tau);
}

class Filter{
public:
    virtual glm::vec2 Radius() const = 0;
    virtual double Evaluate(const glm::vec2& p) const = 0;
    virtual double Integral() const = 0;
};


class BoxFilter : public Filter{
public:
    virtual ~BoxFilter() = default;

    BoxFilter(const glm::vec2& radius = glm::vec2 { 0.5f }) : radius { radius }{}

    glm::vec2 Radius() const final{
        return radius;
    }

    double Evaluate(const glm::vec2& p) const final{
        return std::abs(p.x) <= radius.x && std::abs(p.y) <= radius.y;
    }

    double Integral() const final{
        return 4 * radius.x * radius.y;
    }
private:
    glm::vec2 radius;
};


class GaussianFilter : public Filter{
public:
    virtual ~GaussianFilter() = default;

    GaussianFilter(const glm::vec2& radius = glm::vec2 { 1.5f }, double sigma = 0.5f) : radius { radius }, sigma { sigma }, X { Gaussian<double>(radius.x,sigma) }, Y { Gaussian<double>(radius.y,sigma) }{}

    glm::vec2 Radius() const final{
        return radius;
    }

    double Evaluate(const glm::vec2& p) const final{
        return std::max<double>(0, Gaussian<double>(p.x, sigma) - X) * std::max<double>(0, Gaussian<double>(p.y, sigma) - Y);
    }

    double Integral() const final{
        return (GaussianIntegral<double>(-radius.x, radius.x, sigma) - 2 * radius.x * X) * (GaussianIntegral<double>(-radius.y, radius.y, sigma) - 2 * radius.y * Y);
    }
private:
    glm::vec2 radius;
    double sigma;
    double X;
    double Y;
};

class MitchellFilter : public Filter{
public:
    virtual ~MitchellFilter() = default;

    MitchellFilter(const glm::vec2& radius = glm::vec2 { 1.5f }, double b = 1.0 / 3.0, double c = 1.0 / 3.0) : radius { radius }, b { b }, c { c }{}

    glm::vec2 Radius() const final{
        return radius;
    }

    double Evaluate(const glm::vec2& p) const final{
        return Mitchell(2 * p.x / radius.x) * Mitchell(2 * p.y / radius.y);
    }

    double Integral() const final{
        return radius.x * radius.y / 4.0;
    }
private:
    double Mitchell(double x) const{
        double absX = std::abs(x);
        if(absX <= 1.0){
            return 1.0 / 6.0 * ((12 - 9 * b - 6 * c) * absX * absX * absX + (-18 + 12 * b + 6 * c) * absX * absX + (6 - 2 * b));
        } else if(absX <= 2){
            return 1.0 / 6.0 * ((-b - 6 * c) * absX * absX * absX + (6 * b + 30 * c) * absX * absX + (-12 * b - 48 * c) * absX + (8 * b + 24 * c));
        } else return 0;

    }
    glm::vec2 radius;
    double b;
    double c;
};

class LanczosFilter : public Filter{
public:
    virtual ~LanczosFilter() = default;

    LanczosFilter(const glm::vec2& radius = glm::vec2 { 1.5f }, double tau = 3) : radius(radius), tau(tau){}

    glm::vec2 Radius() const final{
        return radius;
    }
    double Evaluate(const glm::vec2& p) const final{
        return WindowedSinc<double>(p.x, radius.x, tau) * WindowedSinc<double>(p.y, radius.y, tau);
    }

    double Integral() const final{
        double acc = 0;
        int samples = 256 * 256;
        int sqrtSamples = std::sqrt(samples);
        double area = 2 * radius.x * radius.y;
        for(int y = 0;y < sqrtSamples; y++){
            for(int x = 0; x < sqrtSamples; x++){
                glm::vec2 u = glm::vec2 { (x + random_double()), (y + random_double()) } / static_cast<float>(sqrtSamples);
                glm::vec2 p = glm::mix(-radius, radius, u);
                acc += Evaluate(p);
            }
        }
        return area * acc / samples;
    }
private:
    glm::vec2 radius;
    double tau;
};
