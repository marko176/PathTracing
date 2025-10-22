#pragma once
#include "Texture.hpp"
#include "Ray.hpp"
#include "Interaction.hpp"
#include "Random.hpp"
#include <iostream>
#include "Onb.hpp"
#include <numbers>
#include <optional>



template <std::floating_point T>
inline T FresnelDielectric(T cosTheta_i, T eta) {
    cosTheta_i = glm::clamp<T>(cosTheta_i, -1., 1.);

    if (cosTheta_i < 0) {
        eta = 1 / eta;
        cosTheta_i = -cosTheta_i;
    }
 
    float sin2Theta_i = 1 - cosTheta_i*cosTheta_i;
    float sin2Theta_t = sin2Theta_i / (eta*eta);
    if (sin2Theta_t >= 1)
        return 1.f;
    float cosTheta_t = std::sqrt(1 - sin2Theta_t);
    float r_parl = (eta * cosTheta_i - cosTheta_t) / (eta * cosTheta_i + cosTheta_t);
    float r_perp = (cosTheta_i - eta * cosTheta_t) / (cosTheta_i + eta * cosTheta_t);
    return (r_parl*r_parl + r_perp*r_perp) / 2;
}

inline glm::vec3 FresnelSchlick(float cos_theta ,const glm::vec3& F0) {
    return F0 + (glm::vec3(1.0f)-F0) * std::pow(1.0f - cos_theta , 5.0f);
}


enum BxDFFlags : u_int32_t {
    None = 0b0,
    Transmissive = 0b1,
    Specular = 0b10
};

struct BxDFSample {
    glm::vec3 f;
    float pdf;
    u_int32_t flags;

    bool isTransmissive() const {
        return flags & BxDFFlags::Transmissive;
    }

    bool isSpecular() const {
        return flags & BxDFFlags::Specular;
    }
};

class MicrofacetDistribution {
public:
    MicrofacetDistribution(float roughnessX, float roughnessY) : alphaX(roughnessToAlpha(roughnessX)) , alphaY(roughnessToAlpha(roughnessY)) {}

    float lambda(const glm::vec3& w) const {
        float cos2Theta = w.z*w.z;
        if(cos2Theta == 0)return 0;//tan is inf
        float sin2Theta = std::max<float>(0, 1 - cos2Theta);
        float sinTheta = std::sqrt(sin2Theta);
        float cosPhi = sinTheta == 0 ? 1 : glm::clamp(w.x / sinTheta, -1.0f , 1.0f);
        float sinPhi = sinTheta == 0 ? 0 : glm::clamp(w.y / sinTheta, -1.0f , 1.0f);
        float alpha2 = (cosPhi * alphaX)*(cosPhi * alphaX) + (sinPhi * alphaY)*(sinPhi * alphaY);
        return (std::sqrt(1.f + alpha2 * sin2Theta / cos2Theta) - 1.0f )/2.0f;
    }

    float D(const glm::vec3 &wh) const {

        if (wh.z <= 0.0) return 0.0;

        float nx = wh.x;
        float ny = wh.y;
        float nz = wh.z;
        float ax = alphaX;
        float ay = alphaY;

        // denominator base: (nx^2/ax^2 + ny^2/ay^2 + nz^2)
        float inv_ax2 = 1.0 / (ax * ax);
        float inv_ay2 = 1.0 / (ay * ay);

        float base = nx*nx * inv_ax2 + ny*ny * inv_ay2 + nz*nz;
        // protect against extremely small base
        if (base <= 1e-12) return 0.0;

        float denom = std::numbers::pi_v<float> * ax * ay * base * base; // squared
        if (denom <= 0.0) return 0.0;

        return 1.0 / denom;
    }

    float G1(const glm::vec3& w) const {
        return 1 / (1 + lambda(w));
    }

    float G(const glm::vec3& wo, const glm::vec3& wi) const {
        return 1 / (1 + lambda(wo) + lambda(wi));
    }

    bool isSmooth() const {
        return std::max(alphaX,alphaY) < 1e-6;
    }

    glm::vec3 sampleWh(const glm::vec3& wo, const glm::vec2& uv) const {
        bool flip = wo.z < 0;
        glm::vec3 wh = sampleGGXVNDF(flip ? -wo : wo,alphaX,alphaY,uv.x,uv.y);
        if(flip)wh = -wh;
        return wh;
    }

    float PDF(const glm::vec3& wo, const glm::vec3& wh) const {
        return D(wh) * G1(wo) * std::abs(dot(wo, wh)  / wo.z);
    }

    static float roughnessToAlpha(float roughness){
        return roughness*roughness;
    }

protected:
    // Eric Heitz, Sampling the GGX Distribution of Visible Normals, Journal of Computer Graphics Techniques (JCGT), vol. 7, no. 4, 1-13, 2018 (revised 2019-06-17)
    // Available online http://jcgt.org/published/0007/04/01/
    // Input Ve: view direction
    // Input alpha_x, alpha_y: roughness parameters
    // Input U1, U2: uniform random numbers
    // Output Ne: normal sampled with PDF D_Ve(Ne) = G1(Ve) * max(0, dot(Ve, Ne)) * D(Ne) / Ve.z
    glm::vec3 sampleGGXVNDF(const glm::vec3& Ve, float alpha_x, float alpha_y, float U1, float U2) const
    {
        // Section 3.2: transforming the view direction to the hemisphere configuration
        glm::vec3 Vh = glm::normalize(glm::vec3(alpha_x * Ve.x, alpha_y * Ve.y, Ve.z));
        // Section 4.1: orthonormal basis (with special case if cross product is zero)
        float lensq = Vh.x * Vh.x + Vh.y * Vh.y;
        glm::vec3 T1 = lensq > 0 ? glm::vec3(-Vh.y, Vh.x, 0) * glm::inversesqrt(lensq) : glm::vec3(1,0,0);
        glm::vec3 T2 = glm::cross(Vh, T1);
        // Section 4.2: parameterization of the projected area
        float r = std::sqrt(U1);
        float phi = 2.0f * std::numbers::pi_v<float> * U2;
        float t1 = r * std::cos(phi);
        float t2 = r * std::sin(phi);
        float s = 0.5f * (1.0f + Vh.z);
        t2 = (1.0f - s)*std::sqrt(1.0f - t1*t1) + s*t2;
        // Section 4.3: reprojection onto hemisphere
        glm::vec3 Nh = t1*T1 + t2*T2 + std::sqrt(std::max<float>(0.0f, 1.0f - t1*t1 - t2*t2))*Vh;
        // Section 3.4: transforming the normal back to the ellipsoid configuration
        glm::vec3 Ne = glm::normalize(glm::vec3(alpha_x * Nh.x, alpha_y * Nh.y, std::max<float>(0.0f, Nh.z)));
        
        return Ne;
    }
    float alphaX;
    float alphaY;
};




struct Material {
    virtual ~Material() = default;
    virtual std::optional<BxDFSample> scatter(const Ray& ray_in, const SurfaceInteraction& interaction, Ray& scattered,float u,const glm::vec2& UV) const {
        return std::nullopt;
    }

    virtual glm::vec3 emitted(float u,float v) const {
        return glm::vec3(0,0,0);
    }

    virtual bool is_specular(const SurfaceInteraction& interaction) const {
        return false;
    }

    virtual glm::vec3 calc_attenuation(const Ray& r_in, const SurfaceInteraction& interaction, const Ray& scattered) const {
        return {1,1,1};
    }

    virtual float PDF(const Ray& ray_in, const SurfaceInteraction& interaction,const Ray& scattered) const {
        return 0;
    }

    virtual glm::vec3 sample_normalMap(const SurfaceInteraction& interaction) const {
        return interaction.ns;
    }

    virtual std::shared_ptr<Texture> getAlphaMask() const{
        return nullptr;
    }

    virtual float Alpha(float u,float v) const{
        return 1;
    }
};




class lambertian : public Material {
public:
    virtual ~lambertian() = default;
    lambertian(const glm::vec3& albedo) : lambertian(std::make_shared<SolidColor>(albedo)) {}
    lambertian(const std::shared_ptr<Texture>& tex,const std::shared_ptr<Texture>& norm = nullptr,const std::shared_ptr<Texture>& roughnessTexture = std::make_shared<SolidColor>(glm::vec3(1)),const std::shared_ptr<Texture>& metallicTexture = std::make_shared<SolidColor>(glm::vec3(0)),const std::shared_ptr<Texture>& alpha_mask = nullptr) : tex(tex), norm(norm), roughnessTexture(roughnessTexture), metallicTexture(metallicTexture),alpha(alpha_mask) {}

    std::optional<BxDFSample> scatter(const Ray& r_in, const SurfaceInteraction& interaction, Ray& scattered,float u,const glm::vec2& UV) const final {
        //doesnt support smooth material!
        float roughness = GetRoughness(interaction);
        float metallic = metallicTexture->Evaluate(interaction).x;


        onb TBN(glm::dot(r_in.dir,interaction.ns)>0 ? -interaction.ns : interaction.ns); //should be interaction.ns 
        MicrofacetDistribution dist{roughness,roughness};
        float prob = SampleProb(roughness,metallic);
        glm::vec3 wo = TBN.toLocal(-r_in.dir);

        glm::vec3 wi;
        if(u >= prob){
            glm::vec3 wh = dist.sampleWh(wo,UV);
            wi = glm::reflect(-wo,wh);
           
        }else{
            float z = std::sqrt(1.0f - UV.y);

            float phi = 2.0f * std::numbers::pi_v<float> * UV.x;

            float sqrt2 = std::sqrt(UV.y);
            
            float x = std::cos(phi) * sqrt2;
            float y = std::sin(phi) * sqrt2;
            
            wi = {x, y, z};

        }
        if(wi.z<=0){
            //std::cout<<"Not Needed in reflection"<<std::endl;
            //rare case due to rounding off error
            return std::nullopt;
        }
        scattered = Ray(interaction.p, TBN.toWorld(wi),r_in.time);
        
        
        glm::vec3 l = scattered.dir;
        glm::vec3 h = glm::normalize(l - r_in.dir);
        
        float diffuse_pdf = prob * wi.z * std::numbers::inv_pi_v<float>;// max(,0)
        
        glm::vec3 wh = TBN.toLocal(h);
        
        //float D = DD(wh,alpha,alpha);
        
        
        float specular_pdf = (1.0f - prob) * dist.PDF(wo,wh) / (4 * std::abs(glm::dot(wo,wh)));
        float pdf = diffuse_pdf + specular_pdf;
        
        
        
        glm::vec3 textureColor = tex->Evaluate(interaction);
        
        glm::vec3 F0 = glm::mix(glm::vec3(0.04f),textureColor,metallic);
        
        
        float NdotV = glm::clamp(wo.z,0.00001f,1.0f);//wo.z
        
        
        float VdotH = glm::clamp(glm::dot(-r_in.dir, h), 0.0f, 1.0f);
        
        
        
        
        //float G = geometrySmith(NdotV, NdotL, roughness);
        //float G = GG(wo,wi,roughness);
        glm::vec3 F = FresnelSchlick(VdotH, F0);
        
        
        
        glm::vec3 numerator    = dist.D(wh) * dist.G(wo,wi) * F;
        float denominator = 4.0f * NdotV * wi.z;
        
        glm::vec3  specular    = numerator  / denominator;
        
        glm::vec3 kD = (glm::vec3(1.0f) - F) * (1.0f - metallic);
        
        glm::vec3 diffuse = kD * textureColor * std::numbers::inv_pi_v<float>;
        
        glm::vec3 f = diffuse + specular;

        return BxDFSample{f,pdf,BxDFFlags::None};
    }


    float SampleProb(float roughness, float metallic) const{
        if(roughness >= 0.7){
            return 1;
        }else if(metallic <= 0.7 && roughness >= 0.0005f){
            return 0.5;
        }else if(metallic <= 0.9 && roughness >= 0.005f){
            return 0.5;//going higher help but produces fireflys
        }else return roughness * (1.0f - metallic) >= 0.1f ? 0.5f : 0.0f;
    }

    float GetRoughness(const SurfaceInteraction& interaction) const {
        return std::max(roughnessTexture->Evaluate(interaction).x, 0.0005f);
    }

    float PDF(const Ray& ray_in, const SurfaceInteraction& interaction,const Ray& scattered) const final{
        float roughness = GetRoughness(interaction);   
        MicrofacetDistribution dist{roughness,roughness};

        onb TBN(glm::dot(ray_in.dir,interaction.ns)>0 ? -interaction.ns : interaction.ns); 
        glm::vec3 wo = TBN.toLocal(-ray_in.dir);
        glm::vec3 wh = TBN.toLocal(glm::normalize(scattered.dir-ray_in.dir));

        float prob = SampleProb(roughness,metallicTexture->Evaluate(interaction).x);
        
        float diffuse = prob * std::abs(glm::dot(interaction.ns,scattered.dir)) * std::numbers::inv_pi_v<float>;
     
        float specular = dist.PDF(wo,wh) / (4 * std::abs(glm::dot(wo,wh)));

        return diffuse + specular;
    }


    glm::vec3 calc_attenuation(const Ray& r_in, const SurfaceInteraction& interaction, const Ray& scattered) const final {

        glm::vec3 N;
        //<0 == front face
        if(glm::dot(r_in.dir,interaction.ns)>0){
            N = glm::normalize(-interaction.ns);
        }else {
            N = glm::normalize(interaction.ns);
        }
        onb TBN(glm::dot(r_in.dir,interaction.ns)>0 ? -interaction.ns : interaction.ns);  
        glm::vec3 wo = TBN.toLocal(-r_in.dir);
        glm::vec3 wi = TBN.toLocal(scattered.dir);
        float roughness = GetRoughness(interaction);
        MicrofacetDistribution dist{roughness,roughness};                   
        glm::vec3 V = glm::normalize(-r_in.dir);         
        glm::vec3 L = glm::normalize(scattered.dir); 

        

        float metallic = metallicTexture->Evaluate(interaction).x;
        glm::vec3 textureColor = tex->Evaluate(interaction);

        glm::vec3 F0 = glm::mix(glm::vec3(0.04f),textureColor,metallic);
        glm::vec3 H = glm::normalize(V + L);


        float NdotL = glm::dot(N, L);
        float NdotV = glm::clamp(glm::dot(N, V),0.00001f,1.0f);//what if abs?
        //if(NdotL <= 0 )return glm::vec3{0,0,0};

        //float NdotH = glm::clamp(glm::dot(N, H), 0.0f, 1.0f);
        float VdotH = glm::clamp(glm::dot(V, H), 0.0f, 1.0f);

        glm::vec3 F = FresnelSchlick(VdotH, F0);

       


        //onb TBN(interaction);//switch to interaction for multimesh models -> need to put tangent and bitangent on quad intersection!!!!


        glm::vec3 wh = TBN.toLocal(H);

        glm::vec3 numerator = dist.D(wh) * dist.G(wo,wi) * F;
        float denominator = 4.0f * NdotV * NdotL;

        glm::vec3  specular    = numerator  / denominator;
        
        glm::vec3 kD = (glm::vec3(1.0f) - F) * (1.0f - metallic);

        glm::vec3 diffuse = kD * textureColor * std::numbers::inv_pi_v<float>;
    
        return diffuse + specular;
    }

  
    
    std::shared_ptr<Texture> getAlphaMask() const final {
        return alpha;
    }
    float Alpha(float u,float v) const final {
        if(alpha){
            return alpha->Evaluate(SurfaceInteraction({0,0,0},{0,0,0},{u,v})).x;//alpha->color_value(u,v).x
        }
        return tex->alpha(u,v);//should just give 1 channel tex->getChannel( 3 );
    }
    
    glm::vec3 sample_normalMap(const SurfaceInteraction& interaction) const final {
        if(norm == nullptr)return interaction.ns;
        glm::vec3 n_norm = glm::normalize(2.0f * norm->Evaluate(interaction) - glm::vec3(1,1,1));
        return onb(interaction).toWorld(n_norm);
    }

private:
    std::shared_ptr<Texture> norm;
    std::shared_ptr<Texture> alpha;
    std::shared_ptr<Texture> roughnessTexture;
    std::shared_ptr<Texture> metallicTexture;
    std::shared_ptr<Texture> tex;
};

class MicrofacetDielectric : public Material {
public:
    virtual ~MicrofacetDielectric() = default;
    MicrofacetDielectric(float refIndex,const glm::vec3& albedo) : MicrofacetDielectric(refIndex,std::make_shared<SolidColor>(albedo)) {}
    MicrofacetDielectric(float refIndex,float roughness, const glm::vec3& albedo) : MicrofacetDielectric(refIndex,std::make_shared<SolidColor>(albedo),nullptr,std::make_shared<SolidColor>(glm::vec3(roughness))) {}
    MicrofacetDielectric(float refIndex,const std::shared_ptr<Texture>& tex,const std::shared_ptr<Texture>& norm = nullptr,const std::shared_ptr<Texture>& roughnessTexture = std::make_shared<SolidColor>(glm::vec3(0.0)),const std::shared_ptr<Texture>& alpha_mask = nullptr) : ri(refIndex), tex(tex), norm(norm), roughnessTexture(roughnessTexture),alpha(alpha_mask) {}

    std::optional<BxDFSample> scatter(const Ray& r_in, const SurfaceInteraction& interaction, Ray& scattered,float u,const glm::vec2& UV) const final {
        float roughness = GetRoughness(interaction);
        MicrofacetDistribution dist{roughness,roughness};
        onb TBN(interaction);
        glm::vec3 wo = TBN.toLocal(-r_in.dir);
        //this is wrong whould be dot(wo,wh) but it is always > 0
        float eta = glm::dot(-r_in.dir,interaction.ns) > 0 ? 1/ri : ri;//was interaction.ns
        if(dist.isSmooth()){
            glm::vec3 N;
            
            glm::vec3 Ng = glm::dot(r_in.dir,interaction.n)>0 ? -interaction.n : interaction.n;

            //<0 == front face
            if(glm::dot(r_in.dir,interaction.ns)>0){
                N = -interaction.ns;
            }else {
                N = interaction.ns;
            }
            
            float F = FresnelDielectric(glm::dot(-r_in.dir,interaction.ns),ri);//should be r?
    
            float R = F;
            float T = 1.0f - R;
            
            glm::vec3 dir;
            if(u < (R / (R + T))){
                
                dir = TBN.toWorld({-wo.x,-wo.y,wo.z});
               
                glm::vec3 point = r_in.at(interaction.t) + shadowEpsilon * Ng;
                scattered = Ray{point ,dir,r_in.time};
            }else{
                dir = glm::refract(r_in.dir,N,eta);
                
                if(dir == glm::vec3(0,0,0))return std::nullopt;
                glm::vec3 point = r_in.at(interaction.t) - shadowEpsilon * Ng;
                scattered = Ray{point,dir,r_in.time};
            }

            
            return BxDFSample{tex->Evaluate(interaction) / std::abs(glm::dot(interaction.ns,dir)),1,BxDFFlags::Transmissive | BxDFFlags::Specular};
        }

        
        

        

        //bool flip = wo.z < 0;
        //glm::vec3 wh = sampleGGXVNDF(flip ? -wo : wo,alpha,alpha,u1,u2);
        //if(flip)wh = -wh;
        glm::vec3 wh = dist.sampleWh(wo,UV);

        glm::vec3 Ng;//normal on our side


        if(glm::dot(r_in.dir,interaction.n)>0){
            Ng = -interaction.n;
        }else{
            Ng = interaction.n;
        }



        

        //float F = FresnelSchlick(std::abs(glm::dot(wo,wh)),glm::vec3(0.04)).x;
        float F = FresnelDielectric<float>(glm::dot(wo,wh),1/eta);//glm::dot(wo,wh) , wrong becouse it is always > 0
        float R = F;
        float T = 1 - R;
        glm::vec3 wi;

        bool reflect = false;
        if(u < (R / (R + T))){
            reflect = true;
            wi = glm::reflect(-wo, wh);
            if(wo.z * wi.z < 0)return std::nullopt;

            glm::vec3 point = r_in.at(interaction.t) + shadowEpsilon * Ng;
            scattered = Ray{point,TBN.toWorld(wi),r_in.time};
        }else{
            wi = glm::refract(-wo,wh,eta);
            
            if(wo.z * wi.z > 0 || wi.z == 0)return std::nullopt;
            glm::vec3 point = r_in.at(interaction.t) - shadowEpsilon * Ng;
            scattered = Ray{point,TBN.toWorld(wi),r_in.time};
        }  


        float cosTheta_o = wo.z;
        float cosTheta_i = wi.z;
        if (cosTheta_i == 0 || cosTheta_o == 0 ) return std::nullopt;

        if(wh.z<0)wh = -wh;
        
        if (glm::dot(wh, wi) * cosTheta_i < 0.0 || glm::dot(wh, wo) * cosTheta_o < 0.0)
            return std::nullopt;
        
        glm::vec3 textureColor = tex->Evaluate(interaction);

        glm::vec3 f;
        float pdf;
        
        //different results becouse G1(rouhgnes / alpha)
        float p =  dist.PDF(wo,wh);//dist.D(wh) * G1(wo,roughness) / std::abs(wo.z) * std::abs(glm::dot(wo,wh));

        if (reflect) {
            //float p =  dist.D(wh) *G1(wo,roughness) / std::abs(wo.z) * std::abs(glm::dot(wo,wh));
            pdf = p / (4 * std::abs(glm::dot(wo, wh))) * R / (R + T);
            f = textureColor * dist.D(wh) * dist.G(wo,wi) * F / std::abs(4 * cosTheta_i * cosTheta_o);
        } else {    

            //float p =  dist.D(wh) * G1(wo,roughness) / std::abs(wo.z) * std::abs(glm::dot(wo,wh));

            float denom = (glm::dot(wi, wh) + glm::dot(wo, wh) * eta)*(glm::dot(wi, wh) + glm::dot(wo, wh) * eta);
            float dwh_dwi = std::abs(glm::dot(wi, wh)) / denom;

            pdf = p * dwh_dwi /* *T*/ / (R + T);
            float ft =  /*(glm::vec3(1) - F) * */ dist.G(wo,wi) *
                        std::abs(glm::dot(wi, wh) * glm::dot(wo, wh) / (denom * cosTheta_i * cosTheta_o));
            f = textureColor * dist.D(wh) * ft;
        }

        return BxDFSample{f,pdf,BxDFFlags::Transmissive | (roughness < 0.001f ? BxDFFlags::Specular : BxDFFlags::None)};
    }

    float GetRoughness(const SurfaceInteraction& interaction) const {
        return roughnessTexture->Evaluate(interaction).x;
    }

    
    float PDF(const Ray& ray_in, const SurfaceInteraction& interaction,const Ray& scattered) const final{
        float roughness = GetRoughness(interaction);
        MicrofacetDistribution dist{roughness,roughness};
        if(dist.isSmooth())return 0;

        onb TBN(interaction);

        glm::vec3 wo = TBN.toLocal(-ray_in.dir);
        glm::vec3 wi = TBN.toLocal(scattered.dir);

        float cosTheta_o = wo.z;
        float cosTheta_i = wi.z;
        bool reflect = cosTheta_i * cosTheta_o > 0;
        float etap = 1;
        if (!reflect){
            etap = cosTheta_o > 0 ? ri : (1 / ri);
        }
        
        glm::vec3 wh = wi * etap + wo;
        if (cosTheta_i == 0 || cosTheta_o == 0 || glm::dot(wh,wh) == 0) return 0;
        wh = glm::normalize(wh);
        if(wh.z<0)wh = -wh;
        //this causes problems
        if (glm::dot(wh, wi) * cosTheta_i < 0.0 || glm::dot(wh, wo) * cosTheta_o < 0.0)
            return 0;
        
        //float F = FresnelSchlick(std::abs(glm::dot(wo,wh)),glm::vec3(0.04)).x;
        float F = FresnelDielectric<float>(glm::dot(wo,wh),glm::dot(-ray_in.dir,interaction.ns) > 0 ? ri : 1/ri);
        float R = F;
        float T = 1 - R;


        //dist->PDF()
        double pdf = dist.PDF(wo,wh);//DD(wh,alpha,alpha) * G1(wo,roughness) / std::abs(wo.z) * std::abs(glm::dot(wo,wh));
        if (reflect) {
            return pdf / (4 * std::abs(glm::dot(wo, wh))) * R / (R + T);
        } else {    
            double denom = (glm::dot(wi, wh) + glm::dot(wo, wh) / etap)*(glm::dot(wi, wh) + glm::dot(wo, wh) / etap);
            double dwh_dwi = std::abs(glm::dot(wi, wh)) / denom;

            return pdf * dwh_dwi * T / (R + T);
        }
        return 0;
    }

    bool is_specular(const SurfaceInteraction& interaction) const final{
        return true;
    }

    glm::vec3 calc_attenuation(const Ray& r_in, const SurfaceInteraction& interaction, const Ray& scattered) const final {
        
        float roughness = GetRoughness(interaction);
        MicrofacetDistribution dist{roughness,roughness};
        if(dist.isSmooth())return {0,0,0};

        onb TBN(interaction);
        glm::vec3 wo = TBN.toLocal(-r_in.dir);
        glm::vec3 wi = TBN.toLocal(scattered.dir);

        float cosTheta_o = wo.z;
        float cosTheta_i = wi.z;
        bool reflect = cosTheta_i * cosTheta_o > 0;
        float etap = 1;
        if (!reflect){
            etap = cosTheta_o > 0 ? ri : (1 / ri);
        }
        
        glm::vec3 wh = wi * etap + wo;
        if (cosTheta_i == 0 || cosTheta_o == 0 || glm::dot(wh,wh) == 0) return {0,0,0};
        wh = glm::normalize(wh);
        if(wh.z<0)wh = -wh;

        //this causes problems
        if (glm::dot(wh, wi) * cosTheta_i < 0.0 || glm::dot(wh, wo) * cosTheta_o < 0.0)
            return {0,0,0};
       
        float F = FresnelDielectric<float>(glm::dot(wo,wh),glm::dot(-r_in.dir,interaction.ns) > 0 ? ri : 1/ri);
        glm::vec3 textureColor = tex->Evaluate(interaction);
        if (reflect) {
            return textureColor * dist.D(wh) * dist.G(wo,wi) * F / std::abs(4 * cosTheta_i * cosTheta_o);
        } else {    
            double denom = (glm::dot(wi, wh) + glm::dot(wo, wh)/etap) * (glm::dot(wi, wh) + glm::dot(wo, wh)/etap) * cosTheta_i * cosTheta_o;
            float ft =  dist.D(wh) * (1 - F) * dist.G(wo,wi) * std::abs(glm::dot(wi, wh) * glm::dot(wo, wh) / denom);
            return textureColor * ft;
        }
    }

  
    
    std::shared_ptr<Texture> getAlphaMask() const final {
        return alpha;
    }
    float Alpha(float u,float v) const final {
        if(alpha){
            return alpha->Evaluate(SurfaceInteraction({0,0,0},{0,0,0},{u,v})).x;//alpha->color_value(u,v).x
        }
        return tex->alpha(u,v);//should just give 1 channel tex->getChannel( 3 );
    }
    
    glm::vec3 sample_normalMap(const SurfaceInteraction& interaction) const final {
        if(norm == nullptr)return interaction.ns;
        glm::vec3 n_norm = glm::normalize(2.0f * norm->Evaluate(interaction) - glm::vec3(1,1,1));
        return onb(interaction).toWorld(n_norm);
    }
    
private:
    std::shared_ptr<Texture> norm;
    std::shared_ptr<Texture> alpha;
    std::shared_ptr<Texture> roughnessTexture;
    float ri = 1.5;
    std::shared_ptr<Texture> tex;//remove 
};

/*
class dielectric : public Material {
    public:
    dielectric(float r, glm::vec3 color = glm::vec3(1,1,1)) : ri(r) , color(color){}

    

    std::optional<BxDFSample> scatter(const Ray& r_in, const SurfaceInteraction& interaction, Ray& scattered,float u,const glm::vec2& UV) const final{
        float r = ri;
        glm::vec3 N;
        glm::vec3 Ng;//normal on our side
        float cos_theta = glm::dot(r_in.dir,interaction.ns);


        bool entering = glm::dot(r_in.dir,interaction.n)>0;
        if(entering){
            Ng = -interaction.n;
        }else{
            Ng = interaction.n;
            
        }
        //<0 == front face
        if(cos_theta>0){
            N = -interaction.ns;
          
            //Ng = -interaction.n;
            
        }else {
            //front face
            N = interaction.ns;
            //Ng = interaction.n;
            r = 1.0f/r;

        }
        
        float F = FresnelDielectric(glm::dot(-r_in.dir,interaction.ns),ri);//should be r?
 
        float R = F;
        float T = 1.0f - R;
        
        //if tir -> F = 1
        glm::vec3 dir;
        if(u < (R / (R + T))){
       
            dir = glm::reflect(r_in.dir,N);
       
            glm::vec3 point = r_in.at(interaction.t) + shadowEpsilon * Ng;

            scattered = Ray{point ,dir,r_in.time};

        }else{
            dir = glm::refract(r_in.dir,N,r);
            
            if(dir == glm::vec3(0,0,0))return std::nullopt;
            glm::vec3 point = r_in.at(interaction.t) - shadowEpsilon * Ng;
     
            scattered = Ray{point,dir,r_in.time};

        }

        
        return BxDFSample{glm::vec3{1,1,1} / std::abs(glm::dot(interaction.ns,scattered.dir)),1,BxDFFlags::Transmissive | BxDFFlags::Specular};
        
    }

    
    bool is_specular(const SurfaceInteraction& interaction) const final{
        return true;
    }

    virtual glm::vec3 calc_attenuation(const Ray& r_in, const SurfaceInteraction& interaction, const Ray& scattered) const {
        return {0,0,0};
    }

    virtual float PDF(const Ray& ray_in, const SurfaceInteraction& interaction,const Ray& scattered) const {
        return 0;
    }

private:


    glm::vec3 color;
    float ri;
};*/


class SpecularConductor : public Material {
public:
    SpecularConductor(const glm::vec3& albedo) : albedo(albedo) {}
  
    std::optional<BxDFSample> scatter(const Ray& r_in, const SurfaceInteraction& interaction, Ray& scattered,float u,const glm::vec2& UV)const final {
        //fix this needs fresnel term (1-F) 
        scattered = Ray(interaction.p,glm::reflect(r_in.dir,interaction.ns));
        float dot = glm::dot(scattered.dir, interaction.ns);
        if(dot<=0)return std::nullopt;//wont work

        return BxDFSample{FresnelSchlick(glm::dot(interaction.ns,-r_in.dir), albedo) / std::abs(dot),1,BxDFFlags::Specular};
    }

    virtual bool is_specular(const SurfaceInteraction& interaction) const final{
        return true;
    }

    private:
    glm::vec3 albedo;
};

