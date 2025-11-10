#pragma once
#include "Interaction.hpp"
#include "Random.hpp"
#include "Util.hpp"




class Sampler{
public:
    virtual ~Sampler() = default;

    virtual int SamplesPerPixel() const = 0;

    virtual void StartPixelSample(const glm::ivec2& p, int index) = 0;

    virtual double get1D() = 0;

    virtual glm::dvec2 get2D() = 0;

    virtual glm::dvec2 getPixel2D() = 0;

    [[nodiscard]] virtual std::shared_ptr<Sampler> Clone() const = 0;
};

class UniformSampler : public Sampler{
public:
    virtual ~UniformSampler() = default;

    UniformSampler(uint32_t samples) : samples(samples), px(0), py(0){}

    constexpr int SamplesPerPixel() const final{ return samples; }

    void StartPixelSample(const glm::ivec2& p, int index) final{
        px = p.x;
        py = p.y;
    }

    double get1D() final{
        return random_double();
    }

    glm::dvec2 get2D() final{
        return { random_double(),random_double() };
    }

    glm::dvec2 getPixel2D() final{
        return get2D();
    }

    [[nodiscard]] std::shared_ptr<Sampler> Clone() const final{
        return std::make_shared<UniformSampler>(samples);
    }
private:
    uint32_t samples;
    uint32_t px;
    uint32_t py;
};

class StratifiedSampler : public Sampler{
public:
    virtual ~StratifiedSampler() = default;

    StratifiedSampler(uint32_t xSamples, uint32_t ySamples) : xSamples(xSamples), ySamples(ySamples){
        px = 0;
        py = 0;
        sampleIndex = 0;
        dimension = 0;
    }

    constexpr int SamplesPerPixel() const final{ return xSamples * ySamples; }

    void StartPixelSample(const glm::ivec2& p, int index) final{
        px = p.x;
        py = p.y;
        sampleIndex = index;
        dimension = 0;
    }

    double get1D() final{
        uint64_t seed = Hash(px, py, dimension);
        uint64_t stratum = PermutationElement(sampleIndex, SamplesPerPixel(), seed);
        ++dimension;
        return (stratum + random_double()) / (SamplesPerPixel());
    }

    glm::dvec2 get2D() final{
        uint64_t seed = Hash(px, py, dimension);
        uint64_t stratum = PermutationElement(sampleIndex, SamplesPerPixel(), seed);
        dimension += 2;
        int sx = stratum % xSamples;
        int sy = stratum / xSamples;

        // local jitter in cell
        double dx = random_double();
        double dy = random_double();
        return {
            (sx + dx) / double(xSamples),
            (sy + dy) / double(ySamples)
        };
    }

    glm::dvec2 getPixel2D() final{
        return get2D();
    }

    [[nodiscard]] std::shared_ptr<Sampler> Clone() const final{
        return std::make_shared<StratifiedSampler>(xSamples, ySamples);
    }
private:
    uint32_t px;
    uint32_t py;
    uint64_t dimension;
    uint64_t sampleIndex;
    uint16_t xSamples;
    uint16_t ySamples;
};