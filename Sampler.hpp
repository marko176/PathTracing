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

    virtual std::array<glm::vec2,4> get2Dx4f() = 0;

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


    //get1Dx8
    std::array<glm::vec2,4> get2Dx4f() final {
        std::array<float,8> floats = RandomFloatx8();
        std::array<glm::vec2,4> ans;
        for(int i = 0;i<4;i++){
            ans[i]={floats[i*2],floats[i*2+1]};
        }
        return ans;
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
        return (stratum + random_float()) / (SamplesPerPixel());
    }

    glm::dvec2 get2D() final{
        uint64_t seed = Hash(px, py, dimension);
        uint64_t stratum = PermutationElement(sampleIndex, SamplesPerPixel(), seed);
        dimension += 2;
        int sx = stratum % xSamples;
        int sy = stratum / xSamples;

        // local jitter in cell
        double dx = random_float();
        double dy = random_float();
        return {
            (sx + dx) / double(xSamples),
            (sy + dy) / double(ySamples)
        };
    }

    std::array<glm::vec2,4> get2Dx4f() final {

        std::array<float,8> floats = RandomFloatx8();
        std::array<glm::vec2,4> ans;
        for(int i = 0;i<4;i++){
            uint64_t seed = Hash(px, py, dimension);
            uint64_t stratum = PermutationElement(sampleIndex, SamplesPerPixel(), seed);//would be good to move this before randomfloatx8 becouse of the avx slowdown
            
            dimension += 2;

            uint64_t sx = stratum % xSamples;
            uint64_t sy = stratum / xSamples;

            float fx = (sx + floats[i*2]) / static_cast<float>(xSamples);
            float fy = (sy + floats[i*2+1]) / static_cast<float>(ySamples);

            ans[i]=glm::vec2{fx,fy};
        }

        return ans;
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