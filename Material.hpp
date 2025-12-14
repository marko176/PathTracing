#pragma once
#include <iostream>
#include <numbers>
#include <optional>
#include "Texture.hpp"
#include "Ray.hpp"
#include "Interaction.hpp"
#include "Random.hpp"
#include "Onb.hpp"

template <std::floating_point T>
inline T FresnelDielectric(T cosTheta_i, T eta){
    cosTheta_i = glm::clamp<T>(cosTheta_i, -1., 1.);

    if(cosTheta_i < 0){
        eta = 1 / eta;
        cosTheta_i = -cosTheta_i;
    }

    float sin2Theta_i = 1 - cosTheta_i * cosTheta_i;
    float sin2Theta_t = sin2Theta_i / (eta * eta);
    if(sin2Theta_t >= 1)
        return 1.f;
    float cosTheta_t = std::sqrt(1 - sin2Theta_t);
    float r_parl = (eta * cosTheta_i - cosTheta_t) / (eta * cosTheta_i + cosTheta_t);
    float r_perp = (cosTheta_i - eta * cosTheta_t) / (cosTheta_i + eta * cosTheta_t);
    return (r_parl * r_parl + r_perp * r_perp) / 2;
}

inline glm::vec3 FresnelSchlick(float cos_theta, const glm::vec3& F0){
    return F0 + (glm::vec3(1.0f) - F0) * std::pow(1.0f - cos_theta, 5.0f);
}


enum BxDFFlags : uint32_t{
    None = 0b0,
    Transmissive = 0b1,
    Specular = 0b10
};

struct BxDFSample{
    glm::vec3 f;
    float pdf;
    uint32_t flags;//BxDFFlags

    bool isTransmissive() const{
        return flags & BxDFFlags::Transmissive;
    }

    bool isSpecular() const{
        return flags & BxDFFlags::Specular;
    }
};

class MicrofacetDistribution{
public:
    MicrofacetDistribution(float roughnessX, float roughnessY) : alphaX(roughnessToAlpha(roughnessX)), alphaY(roughnessToAlpha(roughnessY)){}

    float lambda(const glm::vec3& w) const{
        float cos2Theta = w.z * w.z;
        if(cos2Theta == 0)return 0;//tan is inf
        float sin2Theta = std::max<float>(0, 1 - cos2Theta);
        float sinTheta = std::sqrt(sin2Theta);
        float cosPhi = sinTheta == 0 ? 1 : glm::clamp(w.x / sinTheta, -1.0f, 1.0f);
        float sinPhi = sinTheta == 0 ? 0 : glm::clamp(w.y / sinTheta, -1.0f, 1.0f);
        float alpha2 = (cosPhi * alphaX) * (cosPhi * alphaX) + (sinPhi * alphaY) * (sinPhi * alphaY);
        return (std::sqrt(1.f + alpha2 * sin2Theta / cos2Theta) - 1.0f) / 2.0f;
    }

    float D(const glm::vec3& wh) const{
        float cos2Theta = wh.z * wh.z;
        if(cos2Theta == 0)return 0;
        float cos4Theta = cos2Theta * cos2Theta;
        float sin2Theta = std::max<float>(0, 1 - cos2Theta);
        float sinTheta = std::sqrt(sin2Theta);
        float cosPhi = sinTheta == 0 ? 1 : glm::clamp(wh.x / sinTheta, -1.0f, 1.0f);
        float sinPhi = sinTheta == 0 ? 0 : glm::clamp(wh.y / sinTheta, -1.0f, 1.0f);
        float e = sin2Theta / cos2Theta * ((cosPhi / alphaX) * (cosPhi / alphaX) + (sinPhi / alphaY) * (sinPhi / alphaY));
        float denom = std::numbers::pi_v<float> *alphaX * alphaY * cos4Theta * (1 + e) * (1 + e);
        if(denom <= 0)return std::numeric_limits<float>::infinity();

        return 1 / denom;
    }

    float G1(const glm::vec3& w) const{
        return 1 / (1 + lambda(w));
    }

    float G(const glm::vec3& wo, const glm::vec3& wi) const{
        return 1 / (1 + lambda(wo) + lambda(wi));
    }

    bool isSmooth() const{
        return std::max(alphaX, alphaY) < 1e-6;
    }

    glm::vec3 sampleWh(const glm::vec3& wo, const glm::vec2& uv) const{
        bool flip = wo.z < 0;
        glm::vec3 wh = sampleGGXVNDF(flip ? -wo : wo, alphaX, alphaY, uv.x, uv.y);
        if(flip)wh = -wh;
        return wh;
    }

    float PDF(const glm::vec3& wo, const glm::vec3& wh) const{
        return D(wh) * G1(wo) * std::abs(dot(wo, wh) / wo.z);
    }

    static float roughnessToAlpha(float roughness){
        return roughness * roughness;
    }

protected:
    // Eric Heitz, Sampling the GGX Distribution of Visible Normals, Journal of Computer Graphics Techniques (JCGT), vol. 7, no. 4, 1-13, 2018 (revised 2019-06-17)
    // Available online http://jcgt.org/published/0007/04/01/
    // Input Ve: view direction
    // Input alpha_x, alpha_y: roughness parameters
    // Input U1, U2: uniform random numbers
    // Output Ne: normal sampled with PDF D_Ve(Ne) = G1(Ve) * max(0, dot(Ve, Ne)) * D(Ne) / Ve.z
    glm::vec3 sampleGGXVNDF(const glm::vec3& Ve, float alpha_x, float alpha_y, float U1, float U2) const{
        // Section 3.2: transforming the view direction to the hemisphere configuration
        glm::vec3 Vh = glm::normalize(glm::vec3(alpha_x * Ve.x, alpha_y * Ve.y, Ve.z));
        // Section 4.1: orthonormal basis (with special case if cross product is zero)
        float lensq = Vh.x * Vh.x + Vh.y * Vh.y;
        glm::vec3 T1 = lensq > 0 ? glm::vec3(-Vh.y, Vh.x, 0) * glm::inversesqrt(lensq) : glm::vec3(1, 0, 0);
        glm::vec3 T2 = glm::cross(Vh, T1);
        // Section 4.2: parameterization of the projected area
        float r = std::sqrt(U1);
        float phi = 2.0f * std::numbers::pi_v<float> *U2;
        float t1 = r * std::cos(phi);
        float t2 = r * std::sin(phi);
        float s = 0.5f * (1.0f + Vh.z);
        t2 = (1.0f - s) * std::sqrt(1.0f - t1 * t1) + s * t2;
        // Section 4.3: reprojection onto hemisphere
        glm::vec3 Nh = t1 * T1 + t2 * T2 + std::sqrt(std::max<float>(0.0f, 1.0f - t1 * t1 - t2 * t2)) * Vh;
        // Section 3.4: transforming the normal back to the ellipsoid configuration
        glm::vec3 Ne = glm::normalize(glm::vec3(alpha_x * Nh.x, alpha_y * Nh.y, std::max<float>(0.0f, Nh.z)));

        return Ne;
    }
    float alphaX;
    float alphaY;
};




struct Material{
    virtual ~Material() = default;
    virtual std::optional<BxDFSample> scatter(const Ray& incoming, const SurfaceInteraction& interaction, Ray& scattered, float u, const glm::vec2& uv) const{
        return std::nullopt;
    }

    virtual glm::vec3 calc_attenuation(const Ray& incoming, const SurfaceInteraction& interaction, const Ray& scattered) const{
        return { 1,1,1 };
    }

    virtual float PDF(const Ray& incoming, const SurfaceInteraction& interaction, const Ray& scattered) const{
        return 0;
    }

    virtual glm::vec3 sample_normalMap(const SurfaceInteraction& interaction) const{
        return interaction.ns;
    }

    virtual bool HasAlpha() const{
        return false;
    }

    //glm::vec2 uv or interaction
    virtual bool Alpha(const glm::vec2& uv) const{
        return true;
    }
};


enum class AlphaMode{
    Opaque,
    Blend,
    Mask
};
struct AlphaTester{
    AlphaTester(AlphaMode mode, float cutoff = 0.5) : mode(mode), cutoff(cutoff){}
    AlphaTester() = default;
    bool operator()(float alpha) const{
        switch(mode){
        case AlphaMode::Opaque:
            return true;
        case AlphaMode::Blend:
            return random_float() < alpha;
        case AlphaMode::Mask:
            return alpha > cutoff;
        default:
            return true;
        }
    }
    AlphaMode mode = AlphaMode::Blend;
    float cutoff = 0.5;
};

class MicrofacetDiffuse : public Material{
public:
    virtual ~MicrofacetDiffuse() = default;
    MicrofacetDiffuse(const glm::vec3& albedo) : MicrofacetDiffuse(std::make_shared<SolidColor>(albedo)){}
    MicrofacetDiffuse(const std::shared_ptr<Texture>& tex, const std::shared_ptr<Texture>& norm = nullptr, const std::shared_ptr<Texture>& roughnessTexture = std::make_shared<SolidColor>(glm::vec3(1)), const std::shared_ptr<Texture>& metallicTexture = std::make_shared<SolidColor>(glm::vec3(0)), const std::shared_ptr<Texture>& alpha_mask = nullptr) : tex(tex), norm(norm), roughnessTexture(roughnessTexture == nullptr ? std::make_shared<SolidColor>(1, 1, 1) : roughnessTexture), metallicTexture(metallicTexture == nullptr ? std::make_shared<SolidColor>(0, 0, 0) : metallicTexture), alpha(alpha_mask){}

    std::optional<BxDFSample> scatter(const Ray& incoming, const SurfaceInteraction& interaction, Ray& scattered, float u, const glm::vec2& uv) const final{
        //doesnt support smooth material!
        float roughness = GetRoughness(interaction);
        

        onb TBN(glm::dot(incoming.dir, interaction.ns) > 0 ? -interaction.ns : interaction.ns);
        MicrofacetDistribution dist { roughness,roughness };
        float prob = SampleProb(roughness);
        glm::vec3 wo = TBN.toLocal(-incoming.dir);
        glm::vec3 wi;
        glm::vec3 wh;
        if(u >= prob){
            wh = dist.sampleWh(wo, uv);
            wi = glm::reflect(-wo, wh);
        } else{
            float z = std::sqrt(1.0f - uv.y);

            float phi = 2.0f * std::numbers::pi_v<float> *uv.x;

            float sqrt2 = std::sqrt(uv.y);

            float x = std::cos(phi) * sqrt2;
            float y = std::sin(phi) * sqrt2;

            wi = { x, y, z };
            wh = glm::normalize(wo + wi);
        }
        if(wi.z <= 0){
            return std::nullopt;
        }
        

        float diffuse_pdf = prob * wi.z * std::numbers::inv_pi_v<float>;

        float specular_pdf = (1.0f - prob) * dist.PDF(wo, wh) / (4 * std::abs(glm::dot(wo, wh)));
        float pdf = diffuse_pdf + specular_pdf;



        glm::vec3 textureColor = tex->Evaluate(interaction);
        float metallic = GetMetallic(interaction);//metallic is in b in gltf
        glm::vec3 F0 = glm::mix(glm::vec3(0.04f), textureColor, metallic);

        glm::vec3 F = FresnelSchlick(glm::dot(wi, wh), F0);


        glm::vec3 numerator = dist.D(wh) * dist.G(wo, wi) * F;
        float denominator = std::abs(4.0f * wo.z * wi.z);
        if(denominator == 0)return std::nullopt;

        glm::vec3  specular = numerator / denominator;

        glm::vec3 kD = (glm::vec3(1.0f) - F) * (1.0f - metallic);

        glm::vec3 diffuse = kD * textureColor * std::numbers::inv_pi_v<float>;

        glm::vec3 f = diffuse + specular;

        scattered = Ray(interaction.p, TBN.toWorld(wi), incoming.time);
        return BxDFSample { f,pdf,BxDFFlags::None };
    }


    static inline float SampleProb(float roughness) {
        return roughness >= 0.7 ? 1 : 0.5;
    }

    float GetRoughness(const SurfaceInteraction& interaction) const{
        return std::max<float>(roughnessTexture->Evaluate(interaction).g,0.0001);//roughness is in g slot
    }

    float GetMetallic(const SurfaceInteraction& interaction) const{
        return metallicTexture->Evaluate(interaction).b;//metallic in b slot
    }

    float PDF(const Ray& incoming, const SurfaceInteraction& interaction, const Ray& scattered) const final{
        float roughness = GetRoughness(interaction);
        MicrofacetDistribution dist { roughness,roughness };

        onb TBN(glm::dot(incoming.dir, interaction.ns) > 0 ? -interaction.ns : interaction.ns);
        glm::vec3 wo = TBN.toLocal(-incoming.dir);
        glm::vec3 wh = TBN.toLocal(glm::normalize(scattered.dir - incoming.dir));

        float prob = SampleProb(roughness);//metallic is in b in gltf

        float diffuse = prob * std::abs(glm::dot(interaction.ns, scattered.dir)) * std::numbers::inv_pi_v<float>;

        float specular = dist.PDF(wo, wh) / (4 * std::abs(glm::dot(wo, wh)));

        return diffuse + specular;
    }


    glm::vec3 calc_attenuation(const Ray& incoming, const SurfaceInteraction& interaction, const Ray& scattered) const final{
        onb TBN(glm::dot(incoming.dir, interaction.ns) > 0 ? -interaction.ns : interaction.ns);
        glm::vec3 wo = TBN.toLocal(-incoming.dir);
        glm::vec3 wi = TBN.toLocal(scattered.dir);
        glm::vec3 wh = glm::normalize(wo + wi);
        float roughness = GetRoughness(interaction);
        float metallic = GetMetallic(interaction);//metallic is in b in gltf

        MicrofacetDistribution dist { roughness,roughness };


        glm::vec3 textureColor = tex->Evaluate(interaction);

        glm::vec3 F0 = glm::mix(glm::vec3(0.04f), textureColor, metallic);

        glm::vec3 F = FresnelSchlick(glm::dot(wi, wh), F0);

        glm::vec3 numerator = dist.D(wh) * dist.G(wo, wi) * F;
        float denominator = std::abs(4.0f * wo.z * wi.z);
        if(denominator == 0)return { 0,0,0 };
        glm::vec3  specular = numerator / denominator;

        glm::vec3 kD = (glm::vec3(1.0f) - F) * (1.0f - metallic);

        glm::vec3 diffuse = kD * textureColor * std::numbers::inv_pi_v<float>;

        return diffuse + specular;
    }



    bool HasAlpha() const final{
        return alphaTester.mode != AlphaMode::Opaque;//must test 
    }

    //tex always has just RGB
    //alpha is handled in geometric primitive?
    bool Alpha(const glm::vec2& uv) const final{
        float a = 0;
        if(alpha){
            a = alpha->Evaluate(SurfaceInteraction({ 0,0,0 }, { 0,0,0 }, uv)).x;//alpha->color_value(u,v).x
        } else a = tex->alpha(uv);//should just give 1 channel tex->getChannel( 3 );
        return alphaTester(a);
    }

    glm::vec3 sample_normalMap(const SurfaceInteraction& interaction) const final{
        if(norm == nullptr)return interaction.ns;
        glm::vec3 n_norm = glm::normalize(2.0f * norm->Evaluate(interaction) - glm::vec3(1, 1, 1));
        return onb(interaction).toWorld(n_norm);
    }

    void setAlphaTester(AlphaTester tester){
        alphaTester = tester;
        if(alpha == nullptr && tex->Channels() != 4)alphaTester.mode = AlphaMode::Opaque;
    }
private:
    std::shared_ptr<Texture> tex;
    std::shared_ptr<Texture> norm;
    std::shared_ptr<Texture> roughnessTexture;
    std::shared_ptr<Texture> metallicTexture;
    std::shared_ptr<Texture> alpha;
    AlphaTester alphaTester;
};

class MicrofacetDielectric : public Material{
public:
    virtual ~MicrofacetDielectric() = default;
    MicrofacetDielectric(float refIndex, const glm::vec3& albedo) : MicrofacetDielectric(refIndex, std::make_shared<SolidColor>(albedo)){}
    MicrofacetDielectric(float refIndex, float roughness, const glm::vec3& albedo) : MicrofacetDielectric(refIndex, std::make_shared<SolidColor>(albedo), nullptr, std::make_shared<SolidColor>(glm::vec3(roughness))){}
    MicrofacetDielectric(float refIndex, const std::shared_ptr<Texture>& tex, const std::shared_ptr<Texture>& norm = nullptr, const std::shared_ptr<Texture>& roughnessTexture = std::make_shared<SolidColor>(glm::vec3(0.0)), const std::shared_ptr<Texture>& alpha_mask = nullptr) : ri(refIndex), tex(tex), norm(norm), roughnessTexture(roughnessTexture != nullptr ? roughnessTexture : std::make_shared<SolidColor>(glm::vec3(0))), alpha(alpha_mask){}

    bool Refract(const glm::vec3& wi, glm::vec3 n, float eta, float* etap, glm::vec3* wt) const{
        float cosTheta_i = glm::dot(n, wi);

        if(cosTheta_i < 0){
            eta = 1 / eta;
            cosTheta_i = -cosTheta_i;
            n = -n;
        }

        float sin2Theta_i = std::max<float>(0, 1 - cosTheta_i * cosTheta_i);
        float sin2Theta_t = sin2Theta_i / (eta * eta);

        float cosTheta_t = std::sqrt(std::max<float>(0, 1 - sin2Theta_t));

        *wt = -wi / eta + (cosTheta_i / eta - cosTheta_t) * glm::vec3(n);

        if(etap)
            *etap = eta;

        return true;
    }

    std::optional<BxDFSample> scatter(const Ray& incoming, const SurfaceInteraction& interaction, Ray& scattered, float u, const glm::vec2& uv) const final{
        float roughness = GetRoughness(interaction);
        MicrofacetDistribution dist { roughness,roughness };
        onb TBN(interaction);
        //eta should be just 1.5
        //onb TBN(glm::dot(incoming.dir,interaction.ns)>0 ? -interaction.ns : interaction.ns);  
        //front_face = glm::dot(incoming.dir,interaction.ns)<0
        glm::vec3 wo = TBN.toLocal(-incoming.dir);
        //this is wrong whould be dot(wo,wh) but it is always > 0
        float eta = glm::dot(-incoming.dir, interaction.ns) > 0 ? 1 / ri : ri;//was interaction.ns
        if(ri == 1 || dist.isSmooth()){
            glm::vec3 N = glm::dot(incoming.dir, interaction.ns) > 0 ? -interaction.ns : interaction.ns;

            glm::vec3 Ng = glm::dot(incoming.dir, interaction.n) > 0 ? -interaction.n : interaction.n;

            float F = FresnelDielectric(wo.z, ri);//should be r?

            float R = F;
            float T = 1.0f - R;

            glm::vec3 dir;
            glm::vec3 f;
            float pdf;
            if(u < (R / (R + T))){

                dir = TBN.toWorld({ -wo.x,-wo.y,wo.z });

                glm::vec3 point = incoming.at(interaction.t) + shadowEpsilon * Ng;
                scattered = Ray { point ,dir,incoming.time };
                f = tex->Evaluate(interaction) * R / std::abs(glm::dot(interaction.ns, dir));
                pdf = R / (R + T);
            } else{

                dir = glm::refract(incoming.dir, N, eta);

                if(dir == glm::vec3(0, 0, 0))return std::nullopt;
                glm::vec3 point = incoming.at(interaction.t) - shadowEpsilon * Ng;
                scattered = Ray { point,dir,incoming.time };
                f = tex->Evaluate(interaction) * T / std::abs(glm::dot(interaction.ns, dir));
                pdf = T / (R + T);
            }


            return BxDFSample { f,pdf,BxDFFlags::Transmissive | BxDFFlags::Specular };
        } else{
            glm::vec3 wh = dist.sampleWh(wo, uv);

            glm::vec3 Ng = glm::dot(incoming.dir, interaction.n) > 0 ? -interaction.n : interaction.n;

            //this here causes energy loss? 1/eta 
            //wh.z > 0 always?
            float F = FresnelDielectric<float>(glm::dot(wo, wh), 1 / eta);//glm::dot(wo,wh) , wrong becouse it is always > 0
            float R = F;
            float T = 1 - R;
            glm::vec3 wi;

            if(u < (R / (R + T))){
                wi = glm::reflect(-wo, wh);
                if(wo.z * wi.z < 0)return std::nullopt;

                glm::vec3 point = incoming.at(interaction.t) + shadowEpsilon * Ng;
                scattered = Ray { point,TBN.toWorld(wi),incoming.time };

                float pdf = dist.PDF(wo, wh) / (4 * std::abs(glm::dot(wo, wh))) * R / (R + T);
                glm::vec3 f = tex->Evaluate(interaction) * dist.D(wh) * dist.G(wo, wi) * R / std::abs(4 * wi.z * wo.z);

                return BxDFSample { f,pdf,BxDFFlags::Transmissive | (roughness < 0.001f ? BxDFFlags::Specular : BxDFFlags::None) };
            } else{
                wi = glm::refract(-wo, wh, eta);//wrong eta?
                //bool tir = !Refract(wo,wh,ri,&eta,&wi);
                if(wo.z * wi.z > 0 || wi.z == 0)return std::nullopt;
                glm::vec3 point = incoming.at(interaction.t) - shadowEpsilon * Ng;
                scattered = Ray { point,TBN.toWorld(wi),incoming.time };

                float denom = (glm::dot(wi, wh) + glm::dot(wo, wh) * eta) * (glm::dot(wi, wh) + glm::dot(wo, wh) * eta);
                float dwh_dwi = std::abs(glm::dot(wi, wh)) / denom;

                float pdf = dist.PDF(wo, wh) * dwh_dwi * T / (R + T);

                float ft = T * dist.D(wh) * dist.G(wo, wi) * std::abs(glm::dot(wi, wh) * glm::dot(wo, wh) / (denom * wi.z * wo.z));
                glm::vec3 f = tex->Evaluate(interaction) * ft;

                return BxDFSample { f,pdf,BxDFFlags::Transmissive | (roughness < 0.001f ? BxDFFlags::Specular : BxDFFlags::None) };
            }
        }
    }

    float GetRoughness(const SurfaceInteraction& interaction) const{
        return roughnessTexture->Evaluate(interaction).y;
    }


    float PDF(const Ray& incoming, const SurfaceInteraction& interaction, const Ray& scattered) const final{
        float roughness = GetRoughness(interaction);
        MicrofacetDistribution dist { roughness,roughness };
        if(ri == 1 || dist.isSmooth())return 0;

        onb TBN(interaction);

        glm::vec3 wo = TBN.toLocal(-incoming.dir);
        glm::vec3 wi = TBN.toLocal(scattered.dir);

        float cosTheta_o = wo.z;
        float cosTheta_i = wi.z;
        bool reflect = cosTheta_i * cosTheta_o > 0;
        float etap = 1;
        if(!reflect){
            etap = cosTheta_o > 0 ? ri : (1 / ri);
        }

        glm::vec3 wh = wi * etap + wo;
        if(glm::dot(wh, wh) == 0) return 0;
        wh = glm::normalize(wh);
        if(wh.z < 0)wh = -wh;

        if(glm::dot(wh, wi) * cosTheta_i <= 0.0 || glm::dot(wh, wo) * cosTheta_o <= 0.0)
            return 0;


        float F = FresnelDielectric<float>(glm::dot(wo, wh), ri);
        float R = F;
        float T = 1 - R;


        //dist->PDF()
        float pdf = dist.PDF(wo, wh);//DD(wh,alpha,alpha) * G1(wo,roughness) / std::abs(wo.z) * std::abs(glm::dot(wo,wh));
        if(reflect){
            return pdf / (4 * std::abs(glm::dot(wo, wh))) * R / (R + T);
        } else{
            float denom = (glm::dot(wi, wh) + glm::dot(wo, wh) / etap) * (glm::dot(wi, wh) + glm::dot(wo, wh) / etap);
            float dwh_dwi = std::abs(glm::dot(wi, wh)) / denom;

            return pdf * dwh_dwi * T / (R + T);
        }
    }

    glm::vec3 calc_attenuation(const Ray& incoming, const SurfaceInteraction& interaction, const Ray& scattered) const final{

        float roughness = GetRoughness(interaction);
        MicrofacetDistribution dist { roughness,roughness };
        if(ri == 1 || dist.isSmooth())return { 0,0,0 };

        onb TBN(interaction);
        glm::vec3 wo = TBN.toLocal(-incoming.dir);
        glm::vec3 wi = TBN.toLocal(scattered.dir);

        float cosTheta_o = wo.z;
        float cosTheta_i = wi.z;
        bool reflect = cosTheta_i * cosTheta_o > 0;
        float etap = 1;
        if(!reflect){
            etap = cosTheta_o > 0 ? ri : (1 / ri);
        }

        glm::vec3 wh = wi * etap + wo;
        if(glm::dot(wh, wh) == 0) return { 0,0,0 };
        wh = glm::normalize(wh);
        if(wh.z < 0)wh = -wh;

        if(glm::dot(wh, wi) * cosTheta_i <= 0.0 || glm::dot(wh, wo) * cosTheta_o <= 0.0)
            return { 0,0,0 };


        float F = FresnelDielectric<float>(glm::dot(wo, wh), ri);
        glm::vec3 textureColor = tex->Evaluate(interaction);
        if(reflect){
            return textureColor * dist.D(wh) * dist.G(wo, wi) * F / std::abs(4 * cosTheta_i * cosTheta_o);
        } else{
            float denom = (glm::dot(wi, wh) + glm::dot(wo, wh) / etap) * (glm::dot(wi, wh) + glm::dot(wo, wh) / etap) * cosTheta_i * cosTheta_o;
            float ft = dist.D(wh) * (1 - F) * dist.G(wo, wi) * std::abs(glm::dot(wi, wh) * glm::dot(wo, wh) / denom);
            return textureColor * ft;
        }
    }



    bool HasAlpha() const final{
        return alphaTester.mode != AlphaMode::Opaque;//must test 
    }

    bool Alpha(const glm::vec2& uv) const final{
        float a = 0;
        if(alpha){
            a = alpha->Evaluate(SurfaceInteraction({ 0,0,0 }, { 0,0,0 }, uv)).x;//alpha->color_value(u,v).x
        } else a = tex->alpha(uv);//should just give 1 channel tex->getChannel( 3 );
        return alphaTester(a);
    }

    glm::vec3 sample_normalMap(const SurfaceInteraction& interaction) const final{
        if(norm == nullptr)return interaction.ns;
        glm::vec3 n_norm = glm::normalize(2.0f * norm->Evaluate(interaction) - glm::vec3(1, 1, 1));
        return onb(interaction).toWorld(n_norm);
    }

    void setAlphaTester(AlphaTester tester){
        alphaTester = tester;
        if(alpha == nullptr && tex->Channels() != 4)alphaTester.mode = AlphaMode::Opaque;
    }

private:
    float ri;
    std::shared_ptr<Texture> tex;
    std::shared_ptr<Texture> norm;
    std::shared_ptr<Texture> roughnessTexture;
    std::shared_ptr<Texture> alpha;
    AlphaTester alphaTester;
};


class ThinDielectric : public Material{
public:
    ThinDielectric(float eta, const std::shared_ptr<Texture>& tex) : ri(eta), albedo(tex != nullptr ? tex : std::make_shared<SolidColor>(glm::vec3(1))){}

    std::optional<BxDFSample> scatter(const Ray& incoming, const SurfaceInteraction& interaction, Ray& scattered, float u, const glm::vec2& uv) const final{
        onb TBN(interaction);
        glm::vec3 wo = TBN.toLocal(-incoming.dir);

        glm::vec3 Ng = glm::dot(incoming.dir, interaction.n) > 0 ? -interaction.n : interaction.n;
        float F = FresnelDielectric(wo.z, ri);

        float R = F;
        float T = 1.0f - R;
        if(R < 1.0f){
            R += T * T * R / (1.0f - R * R);
            T = 1.0f - R;
        }

        glm::vec3 dir;
        glm::vec3 f;
        float pdf;
        if(u < (R / (R + T))){

            dir = TBN.toWorld({ -wo.x,-wo.y,wo.z });

            glm::vec3 point = incoming.at(interaction.t) + shadowEpsilon * Ng;
            scattered = Ray { point ,dir,incoming.time };
            f = glm::vec3(1, 1, 1) * R / std::abs(glm::dot(interaction.ns, dir));
            pdf = R / (R + T);
        } else{

            dir = incoming.dir;

            glm::vec3 point = incoming.at(interaction.t) - shadowEpsilon * Ng;
            scattered = Ray { point,dir,incoming.time };
            f = glm::vec3(1, 1, 1) * T / std::abs(glm::dot(interaction.ns, dir));
            pdf = T / (R + T);
        }

        f *= albedo->Evaluate(interaction);//we need alpha


        return BxDFSample { f,pdf,BxDFFlags::Transmissive | BxDFFlags::Specular };
    }

    glm::vec3 calc_attenuation(const Ray& incoming, const SurfaceInteraction& interaction, const Ray& scattered) const final{
        return { 0,0,0 };
    }

    float PDF(const Ray& incoming, const SurfaceInteraction& interaction, const Ray& scattered) const final{
        return 0;
    }

private:
    float ri;
    std::shared_ptr<Texture> albedo;
};


class SpecularConductor : public Material{
public:
    SpecularConductor(const glm::vec3& albedo) : albedo(albedo){}

    std::optional<BxDFSample> scatter(const Ray& incoming, const SurfaceInteraction& interaction, Ray& scattered, float u, const glm::vec2& uv)const final{
        scattered = Ray(interaction.p, glm::reflect(incoming.dir, interaction.ns), incoming.time);
        float dot = glm::dot(scattered.dir, interaction.ns);
        if(dot <= 0)return std::nullopt;
        return BxDFSample { FresnelSchlick(glm::dot(interaction.ns,-incoming.dir), albedo) / dot,1,BxDFFlags::Specular };
    }

private:
    glm::vec3 albedo;
};

