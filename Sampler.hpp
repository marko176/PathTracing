#pragma once
#include "Hit_record.hpp"
#include "Random.hpp"
#include <cstring>
inline int PermutationElement(uint32_t i, uint32_t l, uint32_t p) {
        uint32_t w = l - 1;
        w |= w >> 1;
        w |= w >> 2;
        w |= w >> 4;
        w |= w >> 8;
        w |= w >> 16;
        do {
            i ^= p;
            i *= 0xe170893d;
            i ^= p >> 16;
            i ^= (i & w) >> 4;
            i ^= p >> 8;
            i *= 0x0929eb3f;
            i ^= p >> 23;
            i ^= (i & w) >> 1;
            i *= 1 | p >> 27;
            i *= 0x6935fa69;
            i ^= (i & w) >> 11;
            i *= 0x74dcb303;
            i ^= (i & w) >> 2;
            i *= 0x9e501cc3;
            i ^= (i & w) >> 2;
            i *= 0xc860a3df;
            i &= w;
            i ^= i >> 5;
        } while (i >= l);
        return (i + p) % l;
    }

inline uint64_t MurmurHash64A(const unsigned char *key, size_t len,
                                           uint64_t seed) {
    const uint64_t m = 0xc6a4a7935bd1e995ull;
    const int r = 47;

    uint64_t h = seed ^ (len * m);

    const unsigned char *end = key + 8 * (len / 8);

    while (key != end) {
        uint64_t k;
        std::memcpy(&k, key, sizeof(uint64_t));
        key += 8;

        k *= m;
        k ^= k >> r;
        k *= m;

        h ^= k;
        h *= m;
    }

    switch (len & 7) {
    case 7:
        h ^= uint64_t(key[6]) << 48;
    case 6:
        h ^= uint64_t(key[5]) << 40;
    case 5:
        h ^= uint64_t(key[4]) << 32;
    case 4:
        h ^= uint64_t(key[3]) << 24;
    case 3:
        h ^= uint64_t(key[2]) << 16;
    case 2:
        h ^= uint64_t(key[1]) << 8;
    case 1:
        h ^= uint64_t(key[0]);
        h *= m;
    };

    h ^= h >> r;
    h *= m;
    h ^= h >> r;

    return h;
}

// Hashing Inline Functions
// http://zimbry.blogspot.ch/2011/09/better-bit-mixing-improving-on.html

inline uint64_t MixBits(uint64_t v) {
    v ^= (v >> 31);
    v *= 0x7fb5d329728ea185;
    v ^= (v >> 27);
    v *= 0x81dadef4bc2dd44d;
    v ^= (v >> 33);
    return v;
}

template <typename T>
inline uint64_t HashBuffer(const T *ptr, size_t size, uint64_t seed = 0) {
    return MurmurHash64A((const unsigned char *)ptr, size, seed);
}

template <typename... Args>
inline uint64_t Hash(Args... args);

template <typename... Args>
inline void hashRecursiveCopy(char *buf, Args...);

template <>
inline void hashRecursiveCopy(char *buf) {}

template <typename T, typename... Args>
inline void hashRecursiveCopy(char *buf, T v, Args... args) {
    memcpy(buf, &v, sizeof(T));
    hashRecursiveCopy(buf + sizeof(T), args...);
}

template <typename... Args>
inline uint64_t Hash(Args... args) {
    // C++, you never cease to amaze: https://stackoverflow.com/a/57246704
    constexpr size_t sz = (sizeof(Args) + ... + 0);
    constexpr size_t n = (sz + 7) / 8;
    uint64_t buf[n];
    hashRecursiveCopy((char *)buf, args...);
    return MurmurHash64A((const unsigned char *)buf, sz, 0);
}

template <typename... Args>
inline float HashFloat(Args... args) {
    return uint32_t(Hash(args...)) * 0x1p-32f;
}




class Sampler {
public:
    virtual int SamplesPerPixel() const = 0;

    virtual void StartPixelSample(const glm::ivec2& p,int index) = 0;

    virtual double get1D() = 0;

    virtual glm::dvec2 get2D() = 0;

    virtual glm::dvec2 GetPixel2D() = 0;

    virtual std::shared_ptr<Sampler> Clone() const = 0;
};

class UniformSampler : public Sampler{
public:
    UniformSampler(uint32_t samples) : px(0), py(0) {

    }

    constexpr int SamplesPerPixel() const final{ return samples; }

    void StartPixelSample(const glm::ivec2& p,int index) final{
        px = p.x;
        py = p.y;
    }

    double get1D() final{
        return random_double();
    }

    glm::dvec2 get2D() final{
        return {random_double(),random_double()};
    }

    glm::dvec2 GetPixel2D() final{
        return get2D();
    }

    std::shared_ptr<Sampler> Clone() const final {
        return std::make_shared<UniformSampler>(samples);
    }
private:
    uint32_t samples;
    uint32_t px;
    uint32_t py;
};

class StratifiedSampler : public Sampler{
public:
    StratifiedSampler(uint32_t xSamples,uint32_t ySamples) : xSamples(xSamples), ySamples(ySamples){
        px = 0;
        py = 0;
        sampleIndex = 0;
        dimension = 0;
    }

    constexpr int SamplesPerPixel() const final { return xSamples*ySamples; }

    void StartPixelSample(const glm::ivec2& p,int index) final {
        px = p.x;
        py = p.y;
        sampleIndex = index;
        dimension = 0;
    }

    double get1D() final {
        uint64_t seed = Hash(px,py, dimension);
        uint64_t stratum = PermutationElement(sampleIndex,SamplesPerPixel(),seed);
        dimension++;
        return (stratum + random_double()) / (SamplesPerPixel());
    }

    glm::dvec2 get2D() final {
        uint64_t seed = Hash(px,py, dimension);
        uint64_t stratum = PermutationElement(sampleIndex,SamplesPerPixel(),seed);
        dimension+=2;
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

    glm::dvec2 GetPixel2D() final {
        return get2D();
    }

    std::shared_ptr<Sampler> Clone() const final {
        return std::make_shared<StratifiedSampler>(xSamples,ySamples);
    }
private:
    uint32_t px;
    uint32_t py;
    uint64_t dimension;
    uint64_t sampleIndex;
    uint16_t xSamples;
    uint16_t ySamples;
};