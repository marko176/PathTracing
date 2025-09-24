#include "Integrators.hpp"
#include "Scene.hpp"
#include "Hit_record.hpp"
#include "Camera.hpp"
#include "Sampler.hpp"
#include <chrono>
#include <thread>

bool Integrator::Unoccluded(const Ray& ray, float t) const {
    return !scene->IntersectPred(ray,t);
}

bool Integrator::Intersect(const Ray& ray, SurfaceInteraction& interaction, float max) const {
    return scene->Intersect(ray,interaction,max);
}

void TileIntegrator::Render() const {
    unsigned int threads = std::thread::hardware_concurrency();
    std::vector<std::thread> workers;
    std::atomic<int> done{0};
    
    constexpr int tileSize = 32;
    glm::ivec2 resolution = camera->GetFilm()->Resolution();
    int tileCount = ((resolution.x + tileSize - 1) / tileSize) * ((resolution.y + tileSize - 1) / tileSize);
    std::mutex consoleMutex;
    uint32_t samples = sampler->SamplesPerPixel();
    auto lamb = [&](){
        int k;
        std::shared_ptr<Sampler> clonedSampler = sampler->Clone();
        while((k = done.fetch_add(1, std::memory_order_relaxed))<tileCount){
            int tileX = k % ((resolution.x + tileSize - 1) / tileSize);//k % 61
            int tileY = k / ((resolution.x + tileSize - 1) / tileSize);
            int minX = tileX * tileSize;
            int minY = tileY * tileSize;
            int maxX = std::min((tileX + 1) * tileSize, resolution.x);
            int maxY = std::min((tileY + 1) * tileSize, resolution.y);
            FilmTile tile = camera->GetFilm()->GetFilmTile({{minX,minY},{maxX,maxY}});
            consoleMutex.lock();
            std::cout<<"\rFinished:"<<std::setw(7)<<std::right<<std::fixed<<std::setprecision(2)<<100 * (done.load())/float(tileCount)<<"%"<<std::flush;
            consoleMutex.unlock();
            for(int y = minY;y < maxY;y++){
                for(int x = minX;x < maxX;x++){
                    VarianceEstimator estimator[3];
                    
                    while(estimator[0].Samples() < 128*samples){//was 32
                        for(int sample_index = 0;sample_index < samples; sample_index++){
                            clonedSampler->StartPixelSample({x,y},sample_index);
                            glm::dvec2 p = glm::dvec2{x,y} + clonedSampler->GetPixel2D();
                            Ray ray = camera->GenerateRay(p,clonedSampler->GetPixel2D(),0);
                            glm::dvec3 color = Li(ray);
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
            camera->GetFilm()->Merge(tile);
        }
    };
    
    auto start = std::chrono::high_resolution_clock::now();
    for(int t = 0;t<threads;t++){
        workers.emplace_back(lamb);
    }
    for(auto& worker : workers)worker.join();
    auto duration = std::chrono::high_resolution_clock::now() - start;
    std::cout<<"\nRender time in ms: "<<std::chrono::duration_cast<std::chrono::milliseconds>(duration).count()<<"\n";
}

glm::vec3 PathIntegrator::Li(Ray ray) const {
    glm::vec3 color = {1,1,1};
    glm::vec3 output = {0,0,0};
    
    int depth = 0;
    int rr_depth = 0;
    float prevPDF = 1;
    bool spec = true;
    while(depth++<maxDepth && (color.x + color.y + color.z) != 0.0f){
        SurfaceInteraction interaction;
        
        if(!Intersect(ray,interaction,1e30f)){
            float a = 0.5f*(ray.dir.y+1.0f);
            return output + color * 1.5f * ((1.0f-a)*glm::vec3(1,0.85,0.55) + a*glm::vec3(0.45,0.65,1));
        }
        
        glm::vec2 random_variables = sampler->get2D();
        glm::vec2 light_random_variables = sampler->get2D();
        float light_selection_random_variable = sampler->get1D();
        float rr_random_variable = sampler->get1D();
        
        glm::vec3 L = {0,0,0};
        if(interaction.AreaLight && (L = interaction.AreaLight->L(interaction,ray)) != glm::vec3(0,0,0)){
            if(spec){
                output += color * L;
            }else{
                float light_pdf = lightSampler->PMF(interaction.AreaLight) * interaction.AreaLight->PDF(interaction,ray);
                float w = prevPDF * prevPDF / (prevPDF * prevPDF + light_pdf * light_pdf);
                output += color * L * w;
            }
        }
        
        Ray new_ray;

        spec = false;
        if(!interaction.mat->scatter(ray,interaction,new_ray,random_variables)){
            return output;//absorbed   
        }
        
        

        if(interaction.mat->is_specular(interaction)){
            color *= interaction.mat->f_PDF(ray,interaction,new_ray);
            ray = new_ray;
            spec = true;
            continue;
        }
        
        float brdfPDF = interaction.mat->PDF(ray,interaction,new_ray);
        prevPDF = brdfPDF;
        output += color * SampleLd(ray,interaction,light_selection_random_variable,light_random_variables);
        glm::vec3 color_attenuation = interaction.mat->calc_attenuation(ray,interaction,new_ray);

        if(brdfPDF <= 0)break;
        color *= color_attenuation / (brdfPDF);
        
        
        if(rr_depth++>3){
            float rr_prob = std::fmin(0.95f,std::fmaxf(std::fmaxf(color.x, color.y), color.z));
            if(rr_random_variable >= rr_prob )break;
            color /= rr_prob;
        }
        ray = new_ray;
        
    }
    return output;
}


glm::vec3 PathIntegrator::SampleLd(const Ray& ray,const SurfaceInteraction& interaction,float u,const glm::vec2& UV) const {
    std::shared_ptr<Light> sampled_light = lightSampler->Sample(u);
    if(sampled_light == nullptr)return {0,0,0};
    LightSample lightSample = sampled_light->sample(UV);
    glm::vec3 lightDir;
    float t = 0;
    if(lightSample.interaction.n == glm::vec3{0,0,0}){
        lightDir = lightSample.dir;
        t = std::numeric_limits<float>::infinity();
    }else{
        lightDir = lightSample.interaction.p - interaction.p;
        t = glm::length(lightDir) - 0.0001f;
    }
    Ray shadow_ray(interaction.p, glm::normalize(lightDir));
    if(glm::dot(interaction.ns,shadow_ray.dir) <= 0 || Unoccluded(shadow_ray,t))return {0,0,0};

    if(sampled_light->isDelta()){
        float light_pdf = lightSampler->PMF(sampled_light) * 1.0f;
        if(light_pdf <= 0)return {0,0,0};
        return sampled_light->L({},shadow_ray) * interaction.mat->calc_attenuation(ray,interaction,shadow_ray) / light_pdf;
    }else{
        float light_pdf = lightSampler->PMF(sampled_light) * sampled_light->PDF(lightSample.interaction,shadow_ray);
        if(light_pdf <= 0)return {0,0,0};
        float w2 = light_pdf*light_pdf;
        float w1 = interaction.mat->PDF(ray,interaction,shadow_ray);
        w1 = w1*w1;
        float w_light = (w2) / (w1 + w2);
        return sampled_light->L(lightSample.interaction,shadow_ray) * interaction.mat->calc_attenuation(ray,interaction,shadow_ray) * w_light / light_pdf;
        
    }
    return {0,0,0};
}

glm::vec3 VolPathIntegrator::Li(Ray ray) const {
    glm::vec3 color = {1,1,1};
    glm::vec3 output = {0,0,0};
    
    int depth = 0;
    int rr_depth = 0;
    float prevPDF = 1;
    bool spec = true;
    while(depth++<maxDepth && (color.x + color.y + color.z) != 0.0f){
        SurfaceInteraction interaction;
        MediumInteraction medInteraction;
        
        if(!Intersect(ray,interaction,std::numeric_limits<float>::infinity())){
            float a = 0.5f*(ray.dir.y+1.0f);
            return output + color * 1.5f * ((1.0f-a)*glm::vec3(1,0.85,0.55) + a*glm::vec3(0.45,0.65,1));
        }

        if(!ray.medium)
            ray.medium = scene->GetMedium();

        if(ray.medium)
            color *= ray.medium->Sample(ray,interaction.t,medInteraction);
        //doesnt work form medium to medium better would be to get medium here?
        
        glm::vec2 random_variables = sampler->get2D();
        glm::vec2 light_random_variables = sampler->get2D();
        float light_selection_random_variable = sampler->get1D();
        float rr_random_variable = sampler->get1D();
        glm::vec2 phase_random_variables = sampler->get2D();
        if(medInteraction.isValid()){
            output += color * SampleLdMedium(ray,medInteraction,light_selection_random_variable,light_random_variables);
            glm::vec3 scattered;
            medInteraction.phaseFunction->Sample(ray.dir,scattered,phase_random_variables);
            ray = Ray(medInteraction.p,scattered,interaction.getMedium(scattered));
            spec = false;
        }else{
            glm::vec3 L = {0,0,0};
            if(interaction.AreaLight && (L = interaction.AreaLight->L(interaction,ray)) != glm::vec3(0,0,0)){
                if(spec){
                    output += color * L;
                }else{
                    float light_pdf = lightSampler->PMF(interaction.AreaLight) * interaction.AreaLight->PDF(interaction,ray);
                    float w = prevPDF * prevPDF / (prevPDF * prevPDF + light_pdf * light_pdf);
                    output += color * L * w;
                }
            }
            spec = false;

            //if mat -> normal
            //if !mat -> fog
            if(!interaction.mat){
                ray.origin = ray.at(interaction.t);
                ray.medium = interaction.getMedium(ray.dir);
                continue;
            }

            Ray new_ray;
    

            
            if(!interaction.mat->scatter(ray,interaction,new_ray,random_variables)){
                return output;//absorbed   
            }

            
            new_ray.medium = interaction.getMedium(new_ray.dir);
    
            if(interaction.mat->is_specular(interaction)){
                color *= interaction.mat->f_PDF(ray,interaction,new_ray);
                ray = new_ray;
                spec = true;
                continue;
            }

            //this helps when medium intersect another object
            //when we bounce we will be in same medium as before
            if(glm::dot(ray.dir,interaction.ns)<=0) 
                new_ray.medium = ray.medium;
            

            float brdfPDF = interaction.mat->PDF(ray,interaction,new_ray);
            prevPDF = brdfPDF;
            output += color * SampleLd(ray,interaction,light_selection_random_variable,light_random_variables);
            glm::vec3 color_attenuation = interaction.mat->calc_attenuation(ray,interaction,new_ray);
    
            
            if(brdfPDF <= 0)break;
            color *= color_attenuation / (brdfPDF);
            ray = new_ray;

            if(rr_depth++>3){
                float rr_prob = std::fmin(0.95f,std::fmaxf(std::fmaxf(color.x, color.y), color.z));

                if(rr_random_variable >= rr_prob)break;
                color /= rr_prob;
            }
        }
        
    }

    return output;
}


glm::vec3 VolPathIntegrator::SampleLd(const Ray& ray,const SurfaceInteraction& interaction,float u,const glm::vec2& UV) const {
    std::shared_ptr<Light> sampled_light = lightSampler->Sample(u);
    if(sampled_light == nullptr)return {0,0,0};
    glm::vec3 Tr = {1,1,1};
    SurfaceInteraction intr;

    LightSample lightSample = sampled_light->sample(UV);
    glm::vec3 lightDir;
    float t = 0;
    if(lightSample.interaction.n == glm::vec3{0,0,0}){
        lightDir = lightSample.dir;
        t = std::numeric_limits<float>::infinity();
    }else{
        lightDir = lightSample.interaction.p - interaction.p;
        t = glm::length(lightDir) - 0.0001f;
    }
    Ray shadow_ray(interaction.p, glm::normalize(lightDir),ray.medium);
    if(glm::dot(interaction.ns,shadow_ray.dir) <= 0 || scene->IntersectTr(shadow_ray,intr,Tr,t))return {0,0,0};
    glm::vec3 f = interaction.mat->calc_attenuation(ray,interaction,shadow_ray);
    if(sampled_light->isDelta()){
        float light_pdf = lightSampler->PMF(sampled_light) * 1.0f;
        if(light_pdf <= 0)return {0,0,0};
        return Tr * sampled_light->L({},shadow_ray) * f / light_pdf;
    }else{
        float light_pdf = lightSampler->PMF(sampled_light) * sampled_light->PDF(lightSample.interaction,shadow_ray);
        if(light_pdf <= 0)return {0,0,0};
        float w2 = light_pdf*light_pdf;
        float w1 = interaction.mat->PDF(ray,interaction,shadow_ray);
        w1 = w1*w1;
        float w_light = (w2) / (w1 + w2);
        return Tr * sampled_light->L(lightSample.interaction,shadow_ray) * f * w_light / light_pdf;
        
    }
    return {0,0,0};
    /*
    if(sampled_light->isDelta()){
        glm::vec3 dir = sampled_light->sample(UV).dir;
        Ray shadow_ray(interaction.p,dir);
        if(glm::dot(interaction.ns,shadow_ray.dir) > 0 && !scene->IntersectTr(shadow_ray,intr, Tr, 1e30f)){
            float light_pdf = lightSampler->PMF(sampled_light);
            if(light_pdf <= 0)return {0,0,0};
            return Tr * sampled_light->L({},shadow_ray) * interaction.mat->calc_attenuation(ray,interaction,shadow_ray) / light_pdf;
        }
    }else{
        LightSample lightSample = sampled_light->sample(UV);
        glm::vec3 to_light = lightSample.interaction.p - interaction.p;
        Ray shadow_ray(interaction.p, glm::normalize(to_light));
        if(glm::dot(interaction.ns,shadow_ray.dir) > 0 && !scene->IntersectTr(shadow_ray,intr, Tr, glm::length(to_light)-0.001f)){
            float light_pdf = lightSampler->PMF(sampled_light) * sampled_light->PDF(lightSample.interaction,shadow_ray);
            if(light_pdf <= 0)return {0,0,0};
            float w2 = light_pdf*light_pdf;
            float w1 = interaction.mat->PDF(ray,interaction,shadow_ray);
            w1 = w1*w1;
            float w_light = (w2) / (w1 + w2);
            return Tr * sampled_light->L(lightSample.interaction,shadow_ray) * interaction.mat->calc_attenuation(ray,interaction,shadow_ray) * w_light / light_pdf;
        }
    }
    return {0,0,0};
    */
}

glm::vec3 VolPathIntegrator::SampleLdMedium(const Ray& ray,const MediumInteraction& interaction,float u,const glm::vec2& UV) const {
    std::shared_ptr<Light> sampled_light = lightSampler->Sample(u);
    if(sampled_light == nullptr)return {0,0,0};
    glm::vec3 Tr = {1,1,1};
    SurfaceInteraction intr;

    LightSample lightSample = sampled_light->sample(UV);
    glm::vec3 lightDir;
    float t = 0;
    if(lightSample.interaction.n == glm::vec3{0,0,0}){
        lightDir = lightSample.dir;
        t = std::numeric_limits<float>::infinity();
    }else{
        lightDir = lightSample.interaction.p - interaction.p;
        t = glm::length(lightDir) - 0.0001f;
    }
    Ray shadow_ray(interaction.p, glm::normalize(lightDir),ray.medium);
    if(scene->IntersectTr(shadow_ray,intr,Tr,t))return {0,0,0};


    float pPhase = interaction.phaseFunction->PDF(ray.dir,shadow_ray.dir);
    glm::vec3 f = glm::vec3(pPhase);
    if(sampled_light->isDelta()){
        float light_pdf = lightSampler->PMF(sampled_light);
        if(light_pdf <= 0)return {0,0,0};
        return Tr * sampled_light->L({},shadow_ray) * f / light_pdf;
    }else{
        float light_pdf = lightSampler->PMF(sampled_light) * sampled_light->PDF(lightSample.interaction,shadow_ray);
        if(light_pdf <= 0)return {0,0,0};
        float w2 = light_pdf*light_pdf;
        float w1 = pPhase;
        w1 = w1*w1;
        float w_light = (w2) / (w1 + w2);
        return Tr * sampled_light->L(lightSample.interaction,shadow_ray) * f * w_light / light_pdf;
        
    }
    return {0,0,0};
}