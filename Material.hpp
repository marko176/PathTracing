#pragma once
#include "Texture.hpp"
#include "Ray.hpp"
#include "Interaction.hpp"
#include "Random.hpp"
#include <iostream>
#include "Onb.hpp"
inline float schlick(float cosine, float ref_idx) {
    float r0 = (1-ref_idx) / (1+ref_idx);
    r0 = r0*r0;
    return r0 + (1-r0)*pow((1 - cosine),5.f);
}


struct Material {
    virtual ~Material() = default;
    virtual bool scatter(const Ray& ray_in, const SurfaceInteraction& interaction, Ray& scattered,const glm::vec2& UV) const {
        return false;
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

    virtual glm::vec3 f_PDF(const Ray& ray_in, const SurfaceInteraction& interaction,Ray& scattered) const {
        return glm::vec3(1,1,1);
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

    static constexpr glm::vec3 fresnel = glm::vec3(0.04);
    glm::vec3 schlick(float cos_theta ,const glm::vec3& F0) const{
        return F0 + (glm::vec3(1)-F0) * std::pow(1.0f - cos_theta , 5.0f);
    }

    float lambda(const glm::vec3& w,float roughness) const {
        if(w.z == 0)return 0;//tan is inf
        return (-1.f + std::sqrt(1.f + roughness*roughness * (w.x*w.x + w.y*w.y)/(w.z*w.z)))/2.0f;
    }

    float G1(const glm::vec3& w, float roughness) const {
        return 1.0f / (1.0f + lambda(w,roughness));
    }


    float distributionGGX(float NdotH,float roughness) const
    {
        float a      = roughness*roughness;
        float a2     = a*a;
        float denom  = (NdotH * NdotH) * (a2 - 1.0f) + 1.0f;
        return a2 * glm::one_over_pi<float>() / (denom * denom);
    }

    float geometrySchlickGGX(float NdotV, float roughness) const
    {
        float r = roughness;
        float k = (r*r) / 2.0f;
        return NdotV / (NdotV * (1.0f - k) + k);
    }

    float geometrySmith(float NdotV, float NdotL, float roughness) const
    {
        float ggx2 = geometrySchlickGGX(NdotV, roughness);
        float ggx1 = geometrySchlickGGX(NdotL, roughness);
        return ggx1 * ggx2;
    }



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

    //std::optional<pdf> -> here pdf is without the D term so no need to 
    bool scatter(const Ray& r_in, const SurfaceInteraction& interaction, Ray& scattered,const glm::vec2& UV) const final {
        
        float u1 = UV.x, u2 = UV.y;
        float roughness = roughnessTexture->Evaluate(interaction).x;
        
        float prob = is_specular(interaction) ? 0 : SampleCosOnly(roughness) ? 1 : 0.5;
        onb TBN(interaction);
        bool sampleVNDF = true;
        if(u1 >= prob){
            u1 = u1*(1.0f + 2.0f*prob) - 2.0f*prob;
            if(!sampleVNDF){
                float cos_theta_h = std::sqrt((1.0f-u1)/(u1*(roughness*roughness - 1.0f) + 1.0f));
                float theta_h = std::acos(cos_theta_h);
                float phi     = 2 * std::numbers::pi_v<float> * u2;
                
                glm::vec3 h_tangent = {
                    sin(theta_h) * cos(phi),
                    sin(theta_h) * sin(phi),
                    cos_theta_h
                };
    
                glm::vec3 h = glm::normalize(TBN.toWorld(h_tangent));
                glm::vec3 l = glm::reflect(r_in.dir, h);
    
                scattered = Ray(interaction.p, glm::normalize(l));

            }else{
                bool flip = glm::dot(interaction.ns,-r_in.dir) < 0;
                glm::vec3 wh = sampleGGXVNDF(flip ? -TBN.toLocal(-r_in.dir) : TBN.toLocal(-r_in.dir),roughness,roughness,u1,u2);
                if(flip)wh = -wh;
                glm::vec3 L = glm::reflect(r_in.dir,TBN.toWorld(wh));

                scattered = Ray(interaction.p, glm::normalize(L));

            }
            
            

        }else{
            /*
            if(prob == 0.5){
                u1 = 2.0f * u1;
            }
            */
            u1 = u1 * (3.0f - 2.f * prob);
            float z = std::sqrt(1.0f - u2);

            float phi = 2.0f * std::numbers::pi_v<float> * u1;

            float sqrt2 = std::sqrt(u2);
            
            float x = std::cos(phi) * sqrt2;
            float y = std::sin(phi) * sqrt2;
            
            glm::vec3 sampled_cosine(x, y, z); // on hemisphere around z+
            //onb uvw(interaction.normal);

            glm::vec3 reflected = TBN.toWorld(sampled_cosine);//was uvw.transform
            scattered = Ray(interaction.p, glm::normalize(reflected));
        }

        if(glm::dot(scattered.dir,interaction.ns)<=0)return false;

        return true;
    }

    bool SampleCosOnly(float roughness) const {
        return roughness > 0.55;
    }


    // anisotropic
    float D_GGX_aniso(const glm::vec3 &h, float ax, float ay) const {
        float x2 = (h.x * h.x)/(ax * ax);
        float y2 = (h.y * h.y)/(ay * ay);
        float d   = x2 + y2 + h.z*h.z;
        return std::numbers::inv_pi_v<float> / (ax * ay * d * d);
    }

    float PDF(const Ray& ray_in, const SurfaceInteraction& interaction,const Ray& scattered) const final{
        glm::vec3 v = -ray_in.dir;
        glm::vec3 l = scattered.dir;
        //glm::vec3 n = rec.normal;
        glm::vec3 h = v+l;
        //float dotVH = glm::dot(v, h);

        
        onb TBN(interaction);
        

        float roughness = roughnessTexture->Evaluate(interaction).x;   
        float prob = is_specular(interaction) ? 0 : SampleCosOnly(roughness) ? 1 : 0.5;
        //float pdf = prob * std::max(glm::dot(interaction.normal,scattered.dir) * std::numbers::inv_pi_v<float>, 0.f) + (1.0f - prob) * distributionGGX(glm::dot(n,h),roughness) * glm::dot(n,h) / (4.0f * std::fabs(dotVH));
        //pdf = prob * std::max(glm::dot(interaction.normal,scattered.dir) * std::numbers::inv_pi_v<float>, 0.f) + (1.0f - prob) * distributionGGX(glm::dot(n,h),roughness) * G1(TBN.toLocal(v),roughness) / (4.0f * std::fabs(dotVH));
        //pdf = prob * std::max(glm::dot(interaction.normal,scattered.dir) * std::numbers::inv_pi_v<float>, 0.f) + (1.0f - prob) * D_GGX_aniso(Ne,roughness*roughness,roughness*roughness) * G1(Ve,roughness) * std::max<float>(0,glm::dot(Ve,Ne))/ Ve.z;
        float diffuse = prob * std::max(glm::dot(interaction.ns,scattered.dir) * std::numbers::inv_pi_v<float>, 0.f);
        float specular = 0;
        if(h != glm::vec3(0) || glm::dot(v,h) == 0){
            float alpha = roughness*roughness;
            glm::vec3 Ve = TBN.toLocal(v);
            glm::vec3 Ne = TBN.toLocal(h);
            h = glm::normalize(h);
            specular = (1.0f - prob) * D_GGX_aniso(Ne,alpha,alpha) * G1(Ve,roughness) * std::max<float>(0,glm::dot(Ve,Ne)) / (4.0f * glm::dot(Ve,Ne) * Ve.z );
        }
        return diffuse + specular;



        //return distributionGGX(glm::dot(n,h),roughness) * glm::dot(n,h) / (4.0f * std::fabs(dotVH));
    }

    

    glm::vec3 f_PDF(const Ray& ray_in, const SurfaceInteraction& interaction,Ray& scattered) const final{
        float pdf = 0;
        glm::vec3 V = -ray_in.dir;
        glm::vec3 L = scattered.dir;
        glm::vec3 N = interaction.ns;
        glm::vec3 H = glm::normalize(V+L);

        
        if (glm::dot(N, L) <= 0.0f || glm::dot(V, H) == 0)
            return glm::vec3(0,0,0);
        onb TBN(interaction);
        glm::vec3 Ve = TBN.toLocal(V);
        glm::vec3 Ne = TBN.toLocal(H);
        float roughness = roughnessTexture->Evaluate(interaction).x;
        pdf = G1(Ve,roughness) * std::max<float>(0,glm::dot(Ve,Ne)) / (4.0f * glm::dot(Ve,Ne) * Ve.z );
        if(pdf == 0 || pdf != pdf)return glm::vec3(0,0,0);



        
        float reflectance = 1;
        //glm::vec3 F0 = fresnel * (reflectance * reflectance);
        float metallic = metallicTexture->Evaluate(interaction).x;
        glm::vec3 textureColor = tex->Evaluate(interaction);
        glm::vec3 F0 = glm::mix(fresnel * (reflectance * reflectance),textureColor,metallic);
        



        float NdotL = glm::clamp(glm::dot(N, L), 0.0f, 1.0f);
        float NdotV = glm::clamp(glm::dot(N, V), 0.0f, 1.0f);
        //float NdotH = glm::clamp(glm::dot(N, H), 0.0f, 1.0f);
        float VdotH = glm::clamp(glm::dot(V, H), 0.0f, 1.0f);

        // 1) D term


        // 2) G term
        float G = geometrySmith(NdotV, NdotL, roughness);

        // 3) F term
        glm::vec3 F = schlick(VdotH, F0);
        
        
        // Specular numerator and denominator
        glm::vec3 numerator    = G * F;
        float      denominator = 4.0f * NdotV * NdotL;
        if(NdotV <= 0){
            return glm::vec3{0,0,0};
        }
        glm::vec3  specular    = numerator / denominator;
        //specular = glm::min(specular,glm::vec3(1000));
        // kS is equal to Fresnel
        glm::vec3 kS = F;
        // kD = diffuse component = 1 − kS, but metallic surfaces have no diffuse
        glm::vec3 kD = (glm::vec3(1.0f) - kS) * (1.0f - metallic);

        // Lambertian diffuse
        glm::vec3 diffuse = kD * textureColor * std::numbers::inv_pi_v<float>;//remove /pi

        // Final BRDF value
        return (diffuse + specular) * NdotL / pdf;//remove ndotl
    }

    bool is_specular(const SurfaceInteraction& interaction) const final {
        float roughness = roughnessTexture->Evaluate(interaction).x;
        return roughness <= 0.02f;
    }
    //f / pdf -> how can we remove D_GGX_aniso(NN,alpha,alpha) from PDF and f ? -> they are always the same ! 
    //

    glm::vec3 calc_attenuation(const Ray& r_in, const SurfaceInteraction& interaction, const Ray& scattered) const final {
        glm::vec3 N = interaction.ns;                               // surface normal ?? maybe geometric normal -> test out!
        glm::vec3 V = glm::normalize(-r_in.dir);          // view vector
        glm::vec3 L = glm::normalize(scattered.dir);        // light/sample vector
        
        float reflectance = 1;
        //glm::vec3 F0 = fresnel * (reflectance * reflectance);
        float metallic = metallicTexture->Evaluate(interaction).x;
        glm::vec3 textureColor = tex->Evaluate(interaction);
        glm::vec3 F0 = glm::mix(fresnel * (reflectance * reflectance),textureColor,metallic);//dont know if needed
        glm::vec3 H = glm::normalize(V + L);

        if(glm::dot(N,L) <= 0 )return glm::vec3{0,0,0};

        float NdotL = glm::clamp(glm::dot(N, L), 0.0f, 1.0f);
        float NdotV = glm::clamp(glm::dot(N, V), 0.0001f, 1.0f);
        //float NdotH = glm::clamp(glm::dot(N, H), 0.0f, 1.0f);
        float VdotH = glm::clamp(glm::dot(V, H), 0.0f, 1.0f);

        // 1) D term
        float roughness = roughnessTexture->Evaluate(interaction).x;
        onb TBN(interaction);//switch to interaction for multimesh models -> need to put tangent and bitangent on quad intersection!!!!


        glm::vec3 NN = TBN.toLocal(H);
        float alpha = roughness*roughness;
        float D = D_GGX_aniso(NN,alpha,alpha);
        //D = distributionGGX(NdotH,roughness);
        // 2) G term
        float G = geometrySmith(NdotV, NdotL, roughness);

        // 3) F term
        glm::vec3 F = schlick(VdotH, F0);
        
        
        // Specular numerator and denominator
        glm::vec3 numerator    = D * G * F;
        float      denominator = 4.0f * NdotV * NdotL;//was max with 1e-7
        //if(NdotV <= 0 || NdotL <= 0){
        //  return glm::vec3{0,0,0};
        //}
        glm::vec3  specular    = numerator  / denominator;
        //specular = glm::min(specular,glm::vec3(1000));
        // kS is equal to Fresnel
        glm::vec3 kS = F;
        // kD = diffuse component = 1 − kS, but metallic surfaces have no diffuse
        glm::vec3 kD = (glm::vec3(1.0f) - kS) * (1.0f - metallic);

        // Lambertian diffuse
        glm::vec3 diffuse = kD * textureColor * std::numbers::inv_pi_v<float>;//remove /pi
    
        // Final BRDF value
        return (diffuse + specular) * NdotL;//remove ndotl
    }

  
    std::shared_ptr<Texture> norm;
    std::shared_ptr<Texture> alpha;
    std::shared_ptr<Texture> roughnessTexture;
    std::shared_ptr<Texture> metallicTexture;
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




        onb TBN(interaction);
        //onb should use tangent bitangent from vertex !
        return TBN.toWorld(n_norm);
        /*
        glm::vec3 N = rec.normal;      // already interpolated & normalized
        glm::vec3 T = rec.tangent;     // already interpolated & normalized
        glm::vec3 B = rec.bitangent;
        // build TBN:
        glm::mat3 TBN = glm::mat3(T, B, N);


        return glm::normalize(TBN * n_norm);*/
    }
    private:
    std::shared_ptr<Texture> tex;
};

/*
class Medium : public Material {
  public:
    Medium(const glm::vec3& albedo, float density = 0.1) : albedo(albedo), d(density) {}

    virtual glm::vec3 calc_attenuation(const Ray& r_in, const SurfaceInteraction& interaction, const Ray& scattered) const override {
        return glm::exp(-interaction.t * (glm::vec3(1) - albedo) * d);//e^0 = 1 -> fully pass
    }

  private:
    glm::vec3 albedo;
    float d;
};
*/


class dielectric : public Material {
    public:
    dielectric(float r, glm::vec3 color = glm::vec3(1,1,1)) : ri(r) , color(color){}
  

    inline glm::vec3 refract(const glm::vec3& vec, const glm::vec3& n, float r) const{
        float cosTheta = std::min(glm::dot(-vec,n),1.0f);
        glm::vec3 outPerp = r * ( vec + cosTheta * n);
        glm::vec3 outParallel = -std::sqrt(std::abs(1.0f - glm::dot(outPerp,outPerp))) * n;
        return outParallel + outPerp;
    }

    bool scatter(const Ray& r_in, const SurfaceInteraction& interaction, Ray& scattered,const glm::vec2& UV) const final{
        /*
        double r = glm::dot(r_in.dir,interaction.ns) < 0 ? (1.0 / ri) : ri;

        double cosTheta = std::min(glm::dot(-r_in.dir,interaction.ns),1.0f);
        double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);

        bool cannotRefract = r * sinTheta > 1.0;
        if(cannotRefract || schlick(cosTheta,r) > UV.x){
            scattered = Ray(interaction.p + 0.002f * N,glm::normalize(glm::reflect(r_in.dir,interaction.ns)));
        }else{
            scattered = Ray(interaction.p - 0.002f * N,glm::normalize(refract(r_in.dir,interaction.ns, (float)r)));
        }

        return true;
        */
        
        float r = ri;
        glm::vec3 N;
        glm::vec3 Ng;//normal on our side
        float cos_theta = glm::dot(r_in.dir,interaction.ns);
        //<0 == front face
        if(cos_theta>0){
            N = -interaction.ns;
            Ng = -interaction.n;
            cos_theta = std::sqrt(1.f - ri*ri*(1.f-cos_theta*cos_theta));
        }else {
            //front face
            N = interaction.ns;
            Ng = interaction.n;
            r = 1.0f/r;
            cos_theta = -cos_theta;
        }
        glm::vec3 dir = glm::refract(r_in.dir,N,r);

        if(dir == glm::vec3(0,0,0) || schlick(cos_theta,r) > UV.x){
            //if reflect and front face

            //case outside to outside -> p += 0
            //case outside to inside -> p += -0.002 n
            //case inside to outside -> p += 0
            //case inside to inside -> p += -0.002 n
            dir = glm::normalize(glm::reflect(r_in.dir,N));
            float eps = glm::dot(interaction.ns,dir) < 0 ? 0.0002f : 0;
            glm::vec3 point = r_in.at(interaction.t) + eps * Ng;
            if(std::abs(glm::dot(interaction.ns,dir))<0.001){
                point = r_in.at(interaction.t) + 0.0002f * dir;//grazing angle
            }
            //should be + in the N direction??
            //also + 0.0005*Ng workes ? but not 0.001
            scattered = Ray{point ,dir};//was + eps*Ng //fix this

        }else{
            
            float eps = false && glm::dot(interaction.ns,dir) < 0 ? 0.0002f : 0;
            glm::vec3 point = r_in.at(interaction.t) - eps * Ng;
            if(std::abs(glm::dot(interaction.ns,dir))<0.001){
                point = r_in.at(interaction.t)+ 0.0002f * dir;//grazing angle
            }
            //case outside to outside -> p += 0 Ng = n
            //case outside to inside -> p += -0.002 Ng = n
            //case inside to outside -> p += -0.002 NG = -n
            //case inside to inside -> p += 0 Ng = -n
            //p is always on ray side (can be backface)
            //refract always into surface
            scattered = Ray{point,dir};//test//was - eps*Ng

        }

        
        return true;
        
    }
    glm::vec3 f_PDF(const Ray& ray_in, const SurfaceInteraction& interaction,Ray& scattered) const final {
        return color;
    }
    
    bool is_specular(const SurfaceInteraction& interaction) const final{
        return true;
    }

    private:

    static double reflectance(double cosine, double refraction_index) {
        // Use Schlick's approximation for reflectance.
        auto r0 = (1 - refraction_index) / (1 + refraction_index);
        r0 = r0*r0;
        return r0 + (1-r0)*std::pow((1 - cosine),5);
    }



    glm::vec3 color;
    float ri;
};

class metal : public Material {
    public:
    metal(const glm::vec3& albedo) : albedo(albedo) {}
  
    bool scatter(const Ray& r_in, const SurfaceInteraction& interaction, Ray& scattered,const glm::vec2& UV)const final {
        scattered = Ray(interaction.p,glm::normalize(glm::reflect(r_in.dir,interaction.ns)));
        return glm::dot(scattered.dir, interaction.ns)>0;//always 1?
    }

    glm::vec3 f_PDF(const Ray& ray_in, const SurfaceInteraction& interaction,Ray& scattered) const final{
        return albedo;
    }
  
    virtual bool is_specular(const SurfaceInteraction& interaction) const final{
        return true;
    }

    virtual glm::vec3 calc_attenuation(const Ray& r_in, const SurfaceInteraction& interaction, const Ray& scattered) const final{
        return albedo;
    }

    private:
    glm::vec3 albedo;
};

