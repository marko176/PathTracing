#include <vector>
#include <iostream>
#include <fstream>
#include <format>
#include <optional>
#include <memory>
#include <stack>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtx/intersect.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <future>
#include <random>
#include "stb_image.h"

#include "Model.hpp"

#include "Texture.hpp"
#include "Material.hpp"
#include "Random.hpp"
#include "Ray.hpp"
#include "Hit_record.hpp"

#include <chrono>


#include "Shape.hpp"
#include "Primitive.hpp"
#include "Light.hpp"
#include "LightSampler.hpp"
#include "Sampler.hpp"
#include "Filter.hpp"
#include "Film.hpp"
#include "Camera.hpp"
#include "Medium.hpp"
#include <OpenImageDenoise/oidn.hpp>
#include "Scene.hpp"
#include "ResourceManager.hpp"

#include "Integrators.hpp"


static constexpr int PrimeTableSize = 1000;
int Primes[PrimeTableSize] = {
    2, 3, 5, 7, 11,
    // Subsequent prime numbers
    13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101,
    103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191,
    193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
    283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389,
    397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491,
    499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607,
    613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719,
    727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829,
    839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953,
    967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051,
    1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153,
    1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259,
    1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367,
    1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471,
    1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567,
    1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 1663,
    1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777,
    1783, 1787, 1789, 1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879,
    1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997, 1999,
    2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063, 2069, 2081, 2083, 2087, 2089, 2099,
    2111, 2113, 2129, 2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, 2221,
    2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, 2293, 2297, 2309, 2311, 2333,
    2339, 2341, 2347, 2351, 2357, 2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417,
    2423, 2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, 2539, 2543, 2549,
    2551, 2557, 2579, 2591, 2593, 2609, 2617, 2621, 2633, 2647, 2657, 2659, 2663, 2671,
    2677, 2683, 2687, 2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741, 2749,
    2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, 2833, 2837, 2843, 2851, 2857,
    2861, 2879, 2887, 2897, 2903, 2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971,
    2999, 3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, 3083, 3089, 3109,
    3119, 3121, 3137, 3163, 3167, 3169, 3181, 3187, 3191, 3203, 3209, 3217, 3221, 3229,
    3251, 3253, 3257, 3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, 3343,
    3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, 3433, 3449, 3457, 3461, 3463,
    3467, 3469, 3491, 3499, 3511, 3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559,
    3571, 3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, 3659, 3671, 3673,
    3677, 3691, 3697, 3701, 3709, 3719, 3727, 3733, 3739, 3761, 3767, 3769, 3779, 3793,
    3797, 3803, 3821, 3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, 3911,
    3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, 4001, 4003, 4007, 4013, 4019,
    4021, 4027, 4049, 4051, 4057, 4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133,
    4139, 4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, 4241, 4243, 4253,
    4259, 4261, 4271, 4273, 4283, 4289, 4297, 4327, 4337, 4339, 4349, 4357, 4363, 4373,
    4391, 4397, 4409, 4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, 4507,
    4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, 4591, 4597, 4603, 4621, 4637,
    4639, 4643, 4649, 4651, 4657, 4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733,
    4751, 4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, 4861, 4871, 4877,
    4889, 4903, 4909, 4919, 4931, 4933, 4937, 4943, 4951, 4957, 4967, 4969, 4973, 4987,
    4993, 4999, 5003, 5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, 5099,
    5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, 5189, 5197, 5209, 5227, 5231,
    5233, 5237, 5261, 5273, 5279, 5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381,
    5387, 5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, 5449, 5471, 5477,
    5479, 5483, 5501, 5503, 5507, 5519, 5521, 5527, 5531, 5557, 5563, 5569, 5573, 5581,
    5591, 5623, 5639, 5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, 5701,
    5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, 5801, 5807, 5813, 5821, 5827,
    5839, 5843, 5849, 5851, 5857, 5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927,
    5939, 5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, 6067, 6073, 6079,
    6089, 6091, 6101, 6113, 6121, 6131, 6133, 6143, 6151, 6163, 6173, 6197, 6199, 6203,
    6211, 6217, 6221, 6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301, 6311,
    6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, 6373, 6379, 6389, 6397, 6421,
    6427, 6449, 6451, 6469, 6473, 6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569,
    6571, 6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, 6679, 6689, 6691,
    6701, 6703, 6709, 6719, 6733, 6737, 6761, 6763, 6779, 6781, 6791, 6793, 6803, 6823,
    6827, 6829, 6833, 6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, 6947,
    6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997, 7001, 7013, 7019, 7027, 7039,
    7043, 7057, 7069, 7079, 7103, 7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193,
    7207, 7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, 7307, 7309, 7321,
    7331, 7333, 7349, 7351, 7369, 7393, 7411, 7417, 7433, 7451, 7457, 7459, 7477, 7481,
    7487, 7489, 7499, 7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, 7573,
    7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, 7649, 7669, 7673, 7681, 7687,
    7691, 7699, 7703, 7717, 7723, 7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823,
    7829, 7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919};








inline uint64_t InverseRadicalInverse(uint64_t inverse, int base,
                                                   int nDigits) {
    uint64_t index = 0;
    for (int i = 0; i < nDigits; ++i) {
        uint64_t digit = inverse % base;
        inverse /= base;
        index = index * base + digit;
    }
    return index;
}

template <typename T>
inline T Mod(T a, T b) {
    T result = a - (a / b) * b;
    return (T)((result < 0) ? result + b : result);
}

inline float RadicalInverse(int baseIndex, uint64_t a) {
    unsigned int base = Primes[baseIndex];
    // We have to stop once reversedDigits is >= limit since otherwise the
    // next digit of |a| may cause reversedDigits to overflow.
    uint64_t limit = ~0ull / base - base;
    float invBase = (float)1 / (float)base, invBaseM = 1;
    uint64_t reversedDigits = 0;
    while (a && reversedDigits < limit) {
        // Extract least significant digit from _a_ and update _reversedDigits_
        uint64_t next = a / base;
        uint64_t digit = a - next * base;
        reversedDigits = reversedDigits * base + digit;
        invBaseM *= invBase;
        a = next;
    }
    return std::min(reversedDigits * invBaseM, 1.0f-std::numeric_limits<float>::epsilon());
}

inline float OwenScrambledRadicalInverse(int baseIndex, uint64_t a,
                                                      uint32_t hash) {
    unsigned int base = Primes[baseIndex];
    // We have to stop once reversedDigits is >= limit since otherwise the
    // next digit of |a| may cause reversedDigits to overflow.
    uint64_t limit = ~0ull / base - base;
    float invBase = (float)1 / (float)base, invBaseM = 1;
    uint64_t reversedDigits = 0;
    int digitIndex = 0;
    while (1 - invBaseM < 1 && reversedDigits < limit) {
        // Compute Owen-scrambled digit for _digitIndex_
        uint64_t next = a / base;
        int digitValue = a - next * base;
        uint32_t digitHash = MixBits(hash ^ reversedDigits);
        digitValue = PermutationElement(digitValue, base, digitHash);
        reversedDigits = reversedDigits * base + digitValue;
        invBaseM *= invBase;
        ++digitIndex;
        a = next;
    }
    return std::min<float>(invBaseM * reversedDigits, 0.99999994f);
}

class HaltonSampler : public Sampler{
  public:
    
    // HaltonSampler Public Methods
    HaltonSampler(int samplesPerPixel, glm::ivec2 fullResolution) : spp(samplesPerPixel) {
        // Find radical inverse base scales and exponents that cover sampling area
        for (int i = 0; i < 2; ++i) {
            int base = (i == 0) ? 2 : 3;
            int scale = 1, exp = 0;
            while (scale < std::min(fullResolution[i], MaxHaltonResolution)) {
                scale *= base;
                ++exp;
            }
            baseScales[i] = scale;
            baseExponents[i] = exp;
        }

        // Compute multiplicative inverses for _baseScales_
        multInverse[0] = multiplicativeInverse(baseScales[1], baseScales[0]);
        multInverse[1] = multiplicativeInverse(baseScales[0], baseScales[1]);
    }

    virtual int SamplesPerPixel() const { return spp; }

    virtual void StartPixelSample(int px,int py, int sampleIndex) {
        //siwcth to ivec
        haltonIndex = 0;
        dimension = 0;
        int sampleStride = baseScales[0] * baseScales[1];
        // Compute Halton sample index for first sample in pixel _p_
        if (sampleStride > 1) {
            glm::ivec2 pm(Mod(px, MaxHaltonResolution), Mod(py, MaxHaltonResolution));
            for (int i = 0; i < 2; ++i) {
                uint64_t dimOffset =
                    (i == 0) ? InverseRadicalInverse(pm[i], 2, baseExponents[i])
                             : InverseRadicalInverse(pm[i], 3, baseExponents[i]);
                haltonIndex +=
                    dimOffset * (sampleStride / baseScales[i]) * multInverse[i];
            }
            haltonIndex %= sampleStride;
        }

        haltonIndex += sampleIndex * sampleStride;
        dimension = 2;
    }

    virtual double get1D() {
        if (dimension >= PrimeTableSize)
            dimension = 2;
        return SampleDimension(dimension++);
    }

    virtual glm::dvec2 get2D() {
        if (dimension + 1 >= PrimeTableSize)
            dimension = 2;
        int dim = dimension;
        dimension += 2;
        return {SampleDimension(dim), SampleDimension(dim + 1)};
    }

    virtual glm::dvec2 GetPixel2D() {
        return {RadicalInverse(0, haltonIndex >> baseExponents[0]),
                RadicalInverse(1, haltonIndex / baseScales[1])};
    }


  private:
    // HaltonSampler Private Methods
    static uint64_t multiplicativeInverse(int64_t a, int64_t n) {
        int64_t x, y;
        extendedGCD(a, n, &x, &y);
        return Mod(x, n);
    }

    static void extendedGCD(uint64_t a, uint64_t b, int64_t *x, int64_t *y) {
        if (b == 0) {
            *x = 1;
            *y = 0;
            return;
        }
        int64_t d = a / b, xp, yp;
        extendedGCD(b, a % b, &xp, &yp);
        *x = yp;
        *y = xp - (d * yp);
    }

    float SampleDimension(int dim) const {
        return OwenScrambledRadicalInverse(dim, haltonIndex,
                                            MixBits(1 + (dim << 4)));
        
    }

    // HaltonSampler Private Members
    int spp;
    static constexpr int MaxHaltonResolution = 128;
    glm::ivec2 baseScales, baseExponents;
    int multInverse[2];
    int64_t haltonIndex = 0;
    int dimension = 0;
};


inline glm::vec3 Li2(Ray curr_ray, const Scene& scene,const std::shared_ptr<Sampler>& sampler,const std::shared_ptr<LightSampler>& LIsampler) {
    glm::vec3 color = {1,1,1};
    glm::vec3 output = {0,0,0};

    
    int depth = 0;
    int rr_depth = 0;
    float prevPDF = 1;
    bool spec = true;

    while(depth++<128 && (color.x + color.y + color.z) != 0.0f){
        SurfaceInteraction interaction;
        MediumInteraction medInteraction;
        
        if(!scene.Intersect(curr_ray,interaction,1e30f)){
            float a = 0.5f*(curr_ray.dir.y+1.0f);
            return output + color * 1.5f * ((1.0f-a)*glm::vec3(1,0.85,0.55) + a*glm::vec3(0.45,0.65,1));
        }

        if(curr_ray.medium){
            color *= curr_ray.medium->Tr(curr_ray,interaction.t);
        }

        
        glm::vec2 random_variables = sampler->get2D();
        glm::vec2 light_random_variables = sampler->get2D();
        float light_selection_random_variable = sampler->get1D();
        float rr_random_variable = sampler->get1D();
        
        glm::vec3 L = {0,0,0};
        if(interaction.AreaLight && (L = interaction.AreaLight->L(interaction,curr_ray)) != glm::vec3(0,0,0)){
            if(spec){
                output += color * L;
            }else{
                double light_pdf = LIsampler->PMF(interaction.AreaLight) * interaction.AreaLight->PDF(interaction,curr_ray);
                float w = prevPDF * prevPDF / (prevPDF * prevPDF + light_pdf * light_pdf);
                output += color * L * w;
            }

        }
        if(interaction.mat){
            Ray new_ray;
    
            //if mat -> normal
            //if !mat -> fog
            
            if(!interaction.mat->scatter(curr_ray,interaction,new_ray,random_variables)){
                return output;//absorbed   
            }
            
            new_ray.medium = interaction.getMedium(new_ray.dir);
    
            if(interaction.mat->is_specular(interaction)){
                color *= interaction.mat->f_PDF(curr_ray,interaction,new_ray);
                curr_ray = new_ray;
                spec = true;
                continue;
            }
            spec = false;

            float brdfPDF = interaction.mat->PDF(curr_ray,interaction,new_ray);
            prevPDF = brdfPDF;
            output += color * LIsampler->SampleLd(curr_ray,interaction,scene.scene_bvh,light_selection_random_variable,light_random_variables);
            glm::vec3 color_attenuation = interaction.mat->calc_attenuation(curr_ray,interaction,new_ray);

            if(brdfPDF <= 0)break;
            color *= color_attenuation / (brdfPDF);
            curr_ray = new_ray;
        }else{
            curr_ray = Ray(curr_ray.at(interaction.t),curr_ray.dir);//we have to have smaller shadow offset for this 0.0001f but for other scenes 0.0005
                                                                    //even 0.0001 is not perfect
            curr_ray.medium = interaction.getMedium(curr_ray.dir);
            spec = false;
        }

        if(rr_depth++>3){
            float rr_prob = std::fmin(0.95f,std::fmaxf(std::fmaxf(color.x, color.y), color.z));
            if(rr_random_variable >= rr_prob )break;
            color /= rr_prob;
        }
        
    }
    return output;
}














void renderPrimFilter(const Scene& scene, int width, int height,std::shared_ptr<LightSampler>& LIsampler, std::string outputImage,const std::shared_ptr<Filter>& filter){

    double fov = 1.7;
    //fov = 0.7;
    fov = 1.7;//dragon


 
    //20 sec
    glm::dvec3 lookfrom = {-1000,300,0};
    lookfrom = {17.3,1.2,7.2};
    //lookfrom = {278,278,-800};
    //lookfrom = {0.3,0.4,1};//dragon
    //lookfrom = {0.3,0.2,1};
    //glm::dvec3 lookfrom = {-900,300,0};
    //glm::dvec3 lookfrom = {1.6,1.6,1.8};
    glm::dvec3 lookat = {0,300,0};
    lookat = {0,0,0};
    double defocus_angle = 0;
    double focus_dist = 6;
    focus_dist = 1;
    //lookat = {278,278,0};
    //lookat = {0,0,0};//dragon
    //glm::dvec3 lookat = {-600,230,-200};

    /*
    double halfWidth  = std::tan(fov * 0.5f);
    double halfHeight = halfWidth * height / width;
    glm::dvec3 up = {0,1,0};
    glm::dvec3 w = glm::normalize(lookfrom-lookat);
    glm::dvec3 u = glm::normalize(glm::cross(up,w));
    glm::dvec3 v = glm::cross(w,u);

    double defocus_radius = focus_dist * std::tan(defocus_angle/2.0f);
    glm::dvec3 defocus_disk_u = u * defocus_radius;
    glm::dvec3 defocus_disk_v = v * defocus_radius;
    */
    unsigned int threads = std::thread::hardware_concurrency();
    std::vector<std::thread> workers;
    std::atomic<int> done{0};
 



    int samples = 25;//64*16*4 -> 4 hours
    int sqrts = std::sqrt(samples);

    std::shared_ptr<Film> film = std::make_shared<Film>(glm::ivec2{width,height},filter);
    constexpr int tileSize = 32;
    int tileCount = ((width + tileSize - 1) / tileSize) * ((height + tileSize - 1) / tileSize);
    std::mutex consoleMutex;
    Camera camera(lookfrom,lookat,fov,film,defocus_angle,focus_dist);
    auto lamb = [&](){
        int k;
        std::shared_ptr<Sampler> sampler = std::make_shared<StratifiedSampler>(sqrts,sqrts);
        while((k = done.fetch_add(1, std::memory_order_relaxed))<tileCount){
            int tileX = k % ((width + tileSize - 1) / tileSize);//k % 61
            int tileY = k / ((width + tileSize - 1) / tileSize);
            int minX = tileX * tileSize;
            int minY = tileY * tileSize;
            int maxX = std::min((tileX + 1) * tileSize, width);
            int maxY = std::min((tileY + 1) * tileSize, height);

            FilmTile tile = camera.GetFilm()->GetFilmTile({{minX,minY},{maxX,maxY}});

            consoleMutex.lock();
            std::cout<<"\rFinished:"<<std::setw(7)<<std::right<<std::fixed<<std::setprecision(2)<<100 * (done.load())/float(tileCount)<<"%"<<std::flush;
            consoleMutex.unlock();

            for(int y = minY;y < maxY;y++){
                for(int x = minX;x < maxX;x++){

                    VarianceEstimator estimator[3];
                    
                    while(estimator[0].Samples() < 128*samples){//was 32
                        for(int sample_index = 0;sample_index < samples; sample_index++){
                            sampler->StartPixelSample({x,y},sample_index);

                            glm::dvec2 p = glm::dvec2{x,y} + sampler->GetPixel2D();
                            Ray ray = camera.GenerateRay(p,sampler->GetPixel2D(),0);

                            glm::dvec3 color = Li2(ray,scene, sampler,LIsampler);
                            if(glm::isnan(color)!=glm::bvec3(false)){
                                std::cout<<"Nan:"<<x<<" "<<y<<"\n";
                                continue;
                            }
                            
                            tile.Add(p,color);
                            color *= glm::dvec3(0.2126f, 0.7152f, 0.0722f);
                            for(int k = 0;k<3;k++)
                                estimator[k].Add(color[k]);
                            
                        }
                        float k = 1.5;//1.5
                        if( estimator[0].RelativeVariance() <= k &&
                            estimator[1].RelativeVariance() <= k &&
                            estimator[2].RelativeVariance() <= k)break;

                    }
            
                }
            }
            camera.GetFilm()->Merge(tile);
        }
    };
    
    auto start = std::chrono::high_resolution_clock::now();
    for(int t = 0;t<threads;t++){
        workers.emplace_back(lamb);
    }
    for(auto& worker : workers)worker.join();

    auto duration = std::chrono::high_resolution_clock::now() - start;
    std::cout<<"\nRender time in ms: "<<std::chrono::duration_cast<std::chrono::milliseconds>(duration).count()<<"\n";
  
    film->WriteImage(outputImage);
}

void temp(){
    auto scene = std::make_shared<Scene>();
    auto white = std::make_shared<Solid_color>(glm::vec3(.9));
    auto green = std::make_shared<Solid_color>(glm::vec3{.2,.3,.1});
    auto light = std::make_shared<lambertian>(glm::vec3(0));
    auto ch =  std::make_shared<lambertian>(glm::vec3{.2,.3,.1});
    auto glass = std::make_shared<dielectric>(1.5,glm::vec3(1));
    ch =  std::make_shared<lambertian>(std::make_shared<CheckerTexture>(white,green,glm::vec2{0.001,0.001}));
    std::shared_ptr<AreaLight> area = std::make_shared<AreaLight>(std::make_shared<QuadShape>(glm::vec3(0.3,1.5,0), glm::vec3(-0.15,0,0), glm::vec3(0,0,-0.15)),glm::vec3(600),false);
    
    scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<QuadShape>(glm::vec3(0.3,1.5,0), glm::vec3(-0.15,0,0), glm::vec3(0,0,-0.15)), light, area, nullptr));//-0.3, -1
    scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<QuadShape>(glm::vec3(-100,-0.3,-100), glm::vec3(1000,0,0), glm::vec3(0,0,1000)), ch, nullptr, nullptr));
    //scene->Add(new Model("/home/markov/Documents/Coding/CPP/raytracing_in_one_weekend/temp_other.assbin"));
    scene->Add(ResourceManager::get_instance().get_model("/home/markov/Documents/Coding/CPP/testing/models/temp_other.assbin",nullptr,nullptr,std::make_shared<HomogeneusMedium>(glm::vec3{0.01f, 0.9f, 0.9f},glm::vec3{1.0f, 0.1f, 0.1f},25.0f)));

    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(0,0.1,-1.2),0.5),std::make_shared<lambertian>(glm::vec3(0.1, 0.2, 0.5)),nullptr));
    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(-1,0,-1),0.5),std::make_shared<dielectric>(1.5),nullptr));
    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(-1,0,-1),0.4),std::make_shared<dielectric>(1/1.5),nullptr));
    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(1,0,-1),0.5),std::make_shared<metal>(glm::vec3(0.8, 0.6, 0.2)),nullptr));
    auto checker = std::make_shared<lambertian>(std::make_shared<CheckerTexture>(white,green,glm::vec2{0.02}));
    scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<SphereShape>(glm::vec3(1,0,-1),0.5),nullptr,nullptr,std::make_shared<HomogeneusMedium>(glm::vec3{0.01f, 0.9f, 0.9f},glm::vec3{1.0f, 0.1f, 0.1f},25.0f)));
        //d_list[1] = new Sphere{glm::vec3(-0.8,1,-0.5), 0.5,
        //                        new Light(glm::vec3(8, 8, 8))};

    std::shared_ptr<LightSampler> ls = std::make_shared<PowerLightSampler>();
    scene->PreProcess();
    ls->Add(scene->GetLights());
    //ls->Add(std::make_shared<PointLight>(glm::vec3(0.3,1.5,0),glm::vec3(6)));
    ls->PreProcess(scene->Bounding_box());
    
    double fov = 1.7;


    glm::dvec3 lookfrom = {-1000,300,0};
    lookfrom = {17.3,1.2,7.2};

    lookfrom = {0.3,0.4,1};//dragon

    glm::dvec3 lookat = {0,300,0};
    lookat = {0,0,0};


    int samples = 64;//64*16*4 -> 4 hours
    int sqrts = std::sqrt(samples);

    std::shared_ptr<Film> film = std::make_shared<Film>(glm::ivec2{1920,1080},std::make_shared<MitchellFilter>());

    auto camera = std::make_shared<Camera>(lookfrom,lookat,fov,film);
    auto sampler = std::make_shared<StratifiedSampler>(sqrts,sqrts);
    auto integrator = std::make_shared<VolPathIntegrator>(scene,camera,sampler,ls,128);

    //renderPrim2(scene,1920,1080,ls,"RenderedScene.ppm");
    integrator->Render();
    camera->GetFilm()->WriteImage("RenderedScene");
   
}

void dragon(){
    Scene scene;
    auto white = std::make_shared<lambertian>(glm::vec3(.73, .73, .73));
    auto light = std::make_shared<lambertian>(glm::vec3(0));
    auto ch =  std::make_shared<lambertian>(glm::vec3{.2,.3,.1});
    
    std::shared_ptr<AreaLight> area = std::make_shared<AreaLight>(std::make_shared<QuadShape>(glm::vec3(0.3,1.5,0), glm::vec3(-0.15,0,0), glm::vec3(0,0,-0.15)),glm::vec3(600),false);
    //what if light color is 0 -> NaN
    scene.Add(std::make_shared<GeometricPrimitive>(std::make_shared<QuadShape>(glm::vec3(0.3,1.5,0), glm::vec3(-0.15,0,0), glm::vec3(0,0,-0.15)), light, area, nullptr));//-0.3, -1
    scene.Add(std::make_shared<GeometricPrimitive>(std::make_shared<QuadShape>(glm::vec3(-100,-0.3,-100), glm::vec3(1000,0,0), glm::vec3(0,0,1000)), ch, nullptr, nullptr));
    //scene.Add(new Model("/home/markov/Documents/Coding/CPP/raytracing_in_one_weekend/temp_other.assbin"));
    scene.Add(ResourceManager::get_instance().get_model("/home/markov/Documents/Coding/CPP/testing/models/temp_other.assbin"));
    
    //scene.Add(new GeometricPrimitive(new SphereShape(glm::vec3(0,0.1,-1.2),0.5),std::make_shared<lambertian>(glm::vec3(0.1, 0.2, 0.5)),nullptr));
    //scene.Add(new GeometricPrimitive(new SphereShape(glm::vec3(-1,0,-1),0.5),std::make_shared<dielectric>(1.5),nullptr));
    //scene.Add(new GeometricPrimitive(new SphereShape(glm::vec3(-1,0,-1),0.4),std::make_shared<dielectric>(1/1.5),nullptr));
    //scene.Add(new GeometricPrimitive(new SphereShape(glm::vec3(1,0,-1),0.5),std::make_shared<metal>(glm::vec3(0.8, 0.6, 0.2)),nullptr));
        //d_list[1] = new Sphere{glm::vec3(-0.8,1,-0.5), 0.5,
        //                        new Light(glm::vec3(8, 8, 8))};

    std::shared_ptr<LightSampler> ls = std::make_shared<PowerLightSampler>();
    scene.PreProcess();
    ls->Add(scene.GetLights());
    ls->PreProcess(scene.Bounding_box());
    

    //renderPrim2(scene,1920,1080,ls,"RenderedScene.ppm");
    renderPrimFilter(scene,1920,1080,ls,"RenderedScene",std::make_shared<MitchellFilter>());
}

int main(){
    
    stbi_set_flip_vertically_on_load(true);
    temp();
    return 0;
    //dragon();
    //return 0;
    Scene scene;
    /*
    auto red   = std::make_shared<lambertian>(glm::vec3(.65, .05, .05));
    auto white = std::make_shared<lambertian>(glm::vec3(.73, .73, .73));
    auto green = std::make_shared<lambertian>(glm::vec3(.12, .45, .15));
    auto light = std::make_shared<lambertian>(glm::vec3(0));
    
    scene.Add(new GeometricPrimitive(new QuadShape(glm::vec3(555,0,0), glm::vec3(0,555,0), glm::vec3(0,0,555)), green, nullptr));
    scene.Add(new GeometricPrimitive(new QuadShape(glm::vec3(0,0,0), glm::vec3(0,555,0), glm::vec3(0,0,555)), red, nullptr));
    std::shared_ptr<AreaLight> area = std::make_shared<AreaLight>(std::make_shared<QuadShape>(glm::vec3(343, 554, 332), glm::vec3(-130,0,0), glm::vec3(0,0,-105)),glm::vec3(15,15,15),false);
    scene.Add(new GeometricPrimitive(new QuadShape(glm::vec3(343, 554, 332), glm::vec3(-130,0,0), glm::vec3(0,0,-105)), light, area));
    scene.Add(new GeometricPrimitive(new QuadShape(glm::vec3(0,0,0), glm::vec3(555,0,0), glm::vec3(0,0,555)), white, nullptr));
    scene.Add(new GeometricPrimitive(new QuadShape(glm::vec3(555,555,555), glm::vec3(-555,0,0), glm::vec3(0,0,-555)), white, nullptr));
    scene.Add(new GeometricPrimitive(new QuadShape(glm::vec3(0,0,555), glm::vec3(555,0,0), glm::vec3(0,555,0)), white, nullptr));
    scene.Add(new GeometricPrimitive(new SphereShape({212,120,147},110),std::make_shared<dielectric>(1.5,glm::vec3(1,1,1)),nullptr,std::make_shared<Medium>(glm::vec3(0.9,0.01,0.01))));
    */
    //std::shared_ptr<LightSampler> ls = std::make_shared<PowerLightSampler>();
    
    
    glm::mat4 pos = glm::mat4(1);
    //pos = glm::translate(pos,{0,-1,0});
    //pos = glm::rotate(pos,glm::radians(-10.0f),glm::vec3(0,1,0));
    //Primitive* modelp = new TransformedPrimitive(std::make_shared<Model>("/home/markov/Documents/Coding/CPP/testing/models/HARD/temp.assbin"),pos);
    scene.Add(ResourceManager::get_instance().get_model("/home/markov/Documents/Coding/CPP/testing/models/HARD/temp.assbin"));
    //scene.Add(new Model("/home/markov/Documents/Coding/CPP/gl/crytek-sponza/temp_other.assbin"));
   
    std::shared_ptr<LightSampler> ls = std::make_shared<PowerLightSampler>();
    ls->Add(std::make_shared<InfiniteLight>(glm::vec3(-1,6,1),25.f*glm::vec3(1,0.93,0.83)));
    /*
    for(int i = 0;i<4;i++){
        for(int j = 0;j<4;j++){
            glm::vec3 color = glm::vec3(3*std::pow((i*4 + j)/15.f,30.f) * 50,3*std::pow((i*4 + j)/15.f,30.f) * 50,3*std::pow((i*4 + j)/15.f,30.f) * 50);
            std::shared_ptr<Shape> shape = std::make_shared<QuadShape>(glm::vec3(7,1.5,1)+glm::vec3(i*2,0,j*2),glm::vec3(0.5,0,0),glm::vec3(0,0,0.5));
            glm::mat4 matrix = glm::mat4(1);
            matrix = glm::translate(matrix,{1,0,0});
            auto lightMaterial = std::make_shared<lambertian>(glm::vec3{0.7});
            std::shared_ptr<AreaLight> area = std::make_shared<AreaLight>(shape,color,true);
            //Primitive* pr = new GeometricPrimitive(new QuadShape(glm::vec3(7,1.5,1)+glm::vec3(i*2,0,j*2),glm::vec3(0.5,0,0),glm::vec3(0,0,0.5)),lightMaterial,area);
            //getLIghts return new Light(transform,light)
            Primitive* pr = new TransformedPrimitive(std::make_shared<GeometricPrimitive>(new QuadShape(glm::vec3(7,1.5,1)+glm::vec3(i*2,0,j*2),glm::vec3(0.5,0,0),glm::vec3(0,0,0.5)),lightMaterial,area),matrix);
            scene.Add(pr);
        }
    }
    */
    /*
    Primitive* pr = new GeometricPrimitive(new SphereShape(glm::vec3(8,0,3),2),
                                            std::make_shared<lambertian>(std::make_shared<Solid_color>(glm::vec3{1,0,0})),nullptr);
    
    scene.Add(pr);
    pr = new GeometricPrimitive(new SphereShape(glm::vec3(8,-202,3),200),
                                            std::make_shared<lambertian>(std::make_shared<Solid_color>(glm::vec3{0,0,1})),nullptr);
    
    scene.Add(pr);
    pr = new GeometricPrimitive(new SphereShape(glm::vec3(7,0,1)+glm::vec3(3*2,0,3*2),1),
                                            std::make_shared<dielectric>(1.5,glm::vec3{1}),nullptr);
    
    scene.Add(pr);
    */
    /*
    Primitive* pr = new GeometricPrimitive(new QuadShape(glm::vec3(7,1,1)+glm::vec3(3*2,0,3*2),glm::vec3(0,1,0),glm::vec3(0,0,1)),
                                            std::make_shared<lambertian>(std::make_shared<Solid_color>(glm::vec3{1}),nullptr,std::make_shared<Solid_color>(glm::vec3{0}),std::make_shared<Solid_color>(glm::vec3{1})),nullptr);
    scene.Add(pr);
    pr = new GeometricPrimitive(new QuadShape(glm::vec3(7,1,1)+glm::vec3(3*2,0,2*2),glm::vec3(0,1,0),glm::vec3(0,0,1)),
                                            std::make_shared<lambertian>(std::make_shared<Solid_color>(glm::vec3{1}),nullptr,std::make_shared<Solid_color>(glm::vec3{0.011}),std::make_shared<Solid_color>(glm::vec3{1})),nullptr);
    scene.Add(pr);   
    pr = new GeometricPrimitive(new QuadShape(glm::vec3(7,1,1)+glm::vec3(3*2,0,1*2),glm::vec3(0,1,0),glm::vec3(0,0,1)),
                                            std::make_shared<lambertian>(std::make_shared<Solid_color>(glm::vec3{1}),nullptr,std::make_shared<Solid_color>(glm::vec3{0.003}),std::make_shared<Solid_color>(glm::vec3{1})),nullptr);
    scene.Add(pr);    
    */
    scene.PreProcess();
    ls->Add(scene.GetLights());
    ls->PreProcess(scene.Bounding_box());
    

    //renderPrim(scene,1920,1080,LIsampler);
    renderPrimFilter(scene,1920,1080,ls,"RenderedScene",std::make_shared<MitchellFilter>());
    //multi_mesh_test_2("/home/markov/Documents/Coding/CPP/testing/models/HARD/temp.assbin");
    //multi_mesh_test_2("/home/markov/Documents/Coding/CPP/gl/crytek-sponza/sponza.obj");


}
