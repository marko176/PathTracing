#include "Integrators.hpp"
#include "Scene.hpp"
#include "Interaction.hpp"
#include "Camera.hpp"
#include "Sampler.hpp"
#include <chrono>
#include <thread>

bool Integrator::Unoccluded(const Ray& ray, float t) const{
    return !scene->IntersectPred(ray, t);
}

bool Integrator::Intersect(const Ray& ray, SurfaceInteraction& interaction, float max) const{
    return scene->Intersect(ray, interaction, max);
}

bool Integrator::IntersectTr(const Ray& ray, SurfaceInteraction& interaction, glm::vec3& Tr, float max) const{
    return scene->IntersectTr(ray, interaction, Tr, max);
}

void TileIntegrator::Render() const{
    unsigned int threads = std::thread::hardware_concurrency();
    std::vector<std::thread> workers;
    std::atomic<int> done { 0 };

    constexpr int tileSize = 32;
    glm::ivec2 resolution = camera->GetFilm()->Resolution();
    int tileCount = ((resolution.x + tileSize - 1) / tileSize) * ((resolution.y + tileSize - 1) / tileSize);
    std::mutex consoleMutex;
    uint32_t samples = sampler->SamplesPerPixel();
    auto lamb = [&](){
        int k;
        std::shared_ptr<Sampler> clonedSampler = sampler->Clone();
        while((k = done.fetch_add(1, std::memory_order_relaxed)) < tileCount){
            int tileX = k % ((resolution.x + tileSize - 1) / tileSize);
            int tileY = k / ((resolution.x + tileSize - 1) / tileSize);
            int minX = tileX * tileSize;
            int minY = tileY * tileSize;
            int maxX = std::min((tileX + 1) * tileSize, resolution.x);
            int maxY = std::min((tileY + 1) * tileSize, resolution.y);
            FilmTile tile = camera->GetFilm()->GetFilmTile({ {minX,minY},{maxX,maxY} });
            consoleMutex.lock();
            std::cout << "\rFinished:" << std::setw(7) << std::right << std::fixed << std::setprecision(2) << 100 * (done.load()) / float(tileCount) << "%" << std::flush;
            consoleMutex.unlock();
            for(int y = minY;y < maxY;y++){
                for(int x = minX;x < maxX;x++){
                    VarianceEstimator estimator[3];

                    while(estimator[0].Samples() < 128 * samples){//was 32
                        for(uint32_t sample_index = 0;sample_index < samples; sample_index++){
                            clonedSampler->StartPixelSample({ x,y }, sample_index);
                            glm::dvec2 p = glm::dvec2 { x,y } + clonedSampler->getPixel2D();
                            float time = clonedSampler->get1D();
                            Ray ray = camera->GenerateRay(p, time, clonedSampler->get2D());
                            glm::dvec3 color = Li(ray);
                            if(glm::isnan(color.x) || glm::isnan(color.y) || glm::isnan(color.z)){
                                std::cout << "Nan:" << x << " " << y << "\n";
                                continue;
                            }

                            tile.Add(p, color);
                            color *= glm::dvec3(0.2126f, 0.7152f, 0.0722f);
                            for(int k = 0;k < 3;k++)
                                estimator[k].Add(color[k]);

                        }
                        float k = 1.5;//lower good for dragon , very low or very high for caustics
                        if(estimator[0].RelativeVariance() <= k &&
                            estimator[1].RelativeVariance() <= k &&
                            estimator[2].RelativeVariance() <= k)break;
                    }

                }
            }
            camera->GetFilm()->Merge(tile);
        }
        };

    auto start = std::chrono::high_resolution_clock::now();
    for(unsigned int t = 0;t < threads;t++){
        workers.emplace_back(lamb);
    }
    for(auto& worker : workers)worker.join();
    auto duration = std::chrono::high_resolution_clock::now() - start;
    std::cout << "\nRender time: " << std::chrono::duration_cast<std::chrono::milliseconds>(duration) << std::endl;
}

glm::vec3 SimplePathIntegrator::Li(Ray ray) const{
    glm::vec3 attenuation = { 1,1,1 };
    glm::vec3 output = { 0,0,0 };

    uint32_t depth = 0;
    uint32_t rr_depth = 0;
    while(depth++ < maxDepth && (attenuation.x + attenuation.y + attenuation.z) != 0.0f){
        SurfaceInteraction interaction;

        if(!Intersect(ray, interaction, std::numeric_limits<float>::infinity())){
            for(auto&& light : scene->infiniteLights){
                output += attenuation * light->Le(ray);
            }
            return output;
        }

        glm::vec2 random_variables = sampler->get2D();
        float scatter_random_variable = sampler->get1D();
        float rr_random_variable = sampler->get1D();

        glm::vec3 L = { 0,0,0 };
        if(interaction.AreaLight && (L = interaction.AreaLight->L(interaction, ray)) != glm::vec3(0, 0, 0)){
            output += attenuation * L;
        }

        Ray new_ray;
        std::optional<BxDFSample> bxdf = interaction.mat->scatter(ray, interaction, new_ray, scatter_random_variable, random_variables);
        if(!bxdf){
            return output;//absorbed   
        }
        new_ray.time = ray.time;

        attenuation *= bxdf->f * std::abs(glm::dot(interaction.ns, new_ray.dir)) / bxdf->pdf;/// brdfPDF;


        if(rr_depth++ > 3){
            float rr_prob = std::fmin(0.95f, std::fmaxf(std::fmaxf(attenuation.x, attenuation.y), attenuation.z));
            if(rr_random_variable >= rr_prob)break;
            attenuation /= rr_prob;
        }
        ray = new_ray;

    }
    return output;
}

glm::vec3 PathIntegrator::Li(Ray ray) const{
    glm::vec3 attenuation = { 1,1,1 };
    glm::vec3 output = { 0,0,0 };

    uint32_t depth = 0;
    uint32_t rr_depth = 0;
    float prevPDF = 1;
    bool spec = true;

    //NEE path splitting count variable -> TODO

    while(depth++ < maxDepth && (attenuation.x + attenuation.y + attenuation.z) != 0.0f){
        SurfaceInteraction interaction;

        if(!Intersect(ray, interaction, std::numeric_limits<float>::infinity())){
            for(auto&& light : scene->infiniteLights){
                glm::vec3 L = light->Le(ray);
                if(spec){
                    output += attenuation * L;
                } else if(prevPDF > 0){
                    float light_pdf = lightSampler->PMF(light) * light->PDF({}, ray);
                    float w = prevPDF * prevPDF / (prevPDF * prevPDF + light_pdf * light_pdf);
                    output += attenuation * L * w;
                }
            }
            return output;
        }

        glm::vec2 random_variables = sampler->get2D();
        glm::vec2 light_random_variables = sampler->get2D();
        float scatter_random_variable = sampler->get1D();
        float light_selection_random_variable = sampler->get1D();
        float rr_random_variable = sampler->get1D();

        glm::vec3 L = { 0,0,0 };
        if(interaction.AreaLight && (L = interaction.AreaLight->L(interaction, ray)) != glm::vec3(0, 0, 0)){
            if(spec){
                output += attenuation * L;
            } else if(prevPDF > 0){
                float light_pdf = lightSampler->PMF(interaction.AreaLight) * interaction.AreaLight->PDF(interaction, ray);
                float w = prevPDF * prevPDF / (prevPDF * prevPDF + light_pdf * light_pdf);
                output += attenuation * L * w;
            }
        }

        if(!interaction.mat){
            spec = true;
            ray.origin = ray.at(interaction.t);
            continue;
        }

        Ray new_ray;



        std::optional<BxDFSample> bxdf = interaction.mat->scatter(ray, interaction, new_ray, scatter_random_variable, random_variables);
        if(!bxdf){
            return output;//absorbed   
        }
        new_ray.time = ray.time;
        if(false && interaction.mat->is_specular(interaction)){

            //is dielectric
            attenuation *= bxdf->f * std::abs(glm::dot(interaction.ns, new_ray.dir)) / bxdf->pdf;

            prevPDF = interaction.mat->PDF(ray, interaction, new_ray);
            spec = false;
            glm::vec3 acc = attenuation * SampleLd(ray, interaction, light_selection_random_variable, light_random_variables);
            output += acc;
            if(acc == glm::vec3(0, 0, 0))spec = true;
            //maybe do same for general reflection not just dielectrics


            ray = new_ray;
            continue;
        }

        spec = bxdf->isSpecular();//is diffuse or is glossy !
        if(!spec){
            output += attenuation * SampleLd(ray, interaction, light_selection_random_variable, light_random_variables);
            prevPDF = interaction.mat->PDF(ray, interaction, new_ray);
        }
        attenuation *= bxdf->f * std::abs(glm::dot(interaction.ns, new_ray.dir)) / bxdf->pdf;

        if(rr_depth++ > 3){
            float rr_prob = std::fmin(0.95f, std::fmaxf(std::fmaxf(attenuation.x, attenuation.y), attenuation.z));
            if(rr_random_variable >= rr_prob)break;
            attenuation /= rr_prob;
        }
        ray = new_ray;

    }
    return output;
}


glm::vec3 PathIntegrator::SampleLd(const Ray& ray, const SurfaceInteraction& interaction, float u, const glm::vec2& UV) const{
    std::shared_ptr<Light> sampled_light = lightSampler->Sample(u);
    if(sampled_light == nullptr)return { 0,0,0 };
    LightSample lightSample = sampled_light->sample(UV, ray.time);
    glm::vec3 lightDir;
    float t;
    if(lightSample.isDeltaInteraction()){
        lightDir = lightSample.dir;
        t = std::numeric_limits<float>::infinity();
    } else{
        lightDir = lightSample.interaction.p - interaction.p;
        t = glm::length(lightDir) - shadowEpsilon;//was 0.0001f
    }
    Ray shadow_ray(interaction.p, glm::normalize(lightDir), ray.time);
    float light_pdf = lightSampler->PMF(sampled_light);

    float dot = glm::dot(interaction.ns, shadow_ray.dir);

    if(light_pdf <= 0 || dot * glm::dot(ray.dir, interaction.ns) >= 0 || !Unoccluded(shadow_ray, t))return { 0,0,0 };

    glm::vec3 f = interaction.mat->calc_attenuation(ray, interaction, shadow_ray) * std::abs(dot);

    if(sampled_light->isDelta()){
        return lightSample.L * f / light_pdf;
    } else{
        light_pdf *= sampled_light->PDF(lightSample.interaction, shadow_ray);
        if(light_pdf <= 0)return { 0,0,0 };
        float w2 = light_pdf * light_pdf;
        float w1 = interaction.mat->PDF(ray, interaction, shadow_ray);
        w1 = w1 * w1;
        float w_light = (w2) / (w1 + w2);
        return sampled_light->L(lightSample.interaction, shadow_ray) * f * w_light / light_pdf;

    }
    return { 0,0,0 };
}

glm::vec3 VolPathIntegrator::Li(Ray ray) const{
    glm::vec3 attenuation = { 1,1,1 };
    glm::vec3 output = { 0,0,0 };

    uint32_t depth = 0;
    uint32_t rr_depth = 0;
    float prevPDF = 1;
    bool spec = true;
    bool EnableCaustics = false;
    while(depth++ < maxDepth && (attenuation.x + attenuation.y + attenuation.z) != 0.0f){
        SurfaceInteraction interaction;
        MediumInteraction medInteraction;

        if(!Intersect(ray, interaction, std::numeric_limits<float>::infinity())){
            for(auto&& light : scene->infiniteLights){
                glm::vec3 L = light->Le(ray);
                if(spec){
                    output += attenuation * L;
                } else if(prevPDF > 0){
                    float light_pdf = lightSampler->PMF(light) * light->PDF({}, ray);
                    float w = prevPDF * prevPDF / (prevPDF * prevPDF + light_pdf * light_pdf);
                    output += attenuation * L * w;
                }
            }
            return output;
        }

        //maybe set medium here
        //


        //test if this is correct?
        //it doesnt work if medium clips another object -> if we are in medium we see we hit outside -> no medium
        //ray.medium = interaction.getInverseMedium(ray.dir);

        if(!ray.medium)
            ray.medium = scene->GetMedium();

        if(ray.medium)
            attenuation *= ray.medium->Sample(ray, interaction.t, medInteraction);
        //doesnt work form medium to medium better would be to get medium here?

        glm::vec2 random_variables = sampler->get2D();
        glm::vec2 light_random_variables = sampler->get2D();
        float scatter_random_variable = sampler->get1D();
        float light_selection_random_variable = sampler->get1D();
        float rr_random_variable = sampler->get1D();
        glm::vec2 phase_random_variables = sampler->get2D();

        //tjedan prije u ponidiljak

        if(medInteraction.isValid()){
            output += attenuation * SampleLd(ray, medInteraction, light_selection_random_variable, light_random_variables);
            output += attenuation * ray.medium->Le();
            glm::vec3 scattered;
            medInteraction.phaseFunction->Sample(ray.dir, scattered, phase_random_variables);
            ray = Ray(medInteraction.p, scattered, ray.time, interaction.getMedium(scattered));
            spec = false;
        } else{
            glm::vec3 L = { 0,0,0 };
            if(interaction.AreaLight && (L = interaction.AreaLight->L(interaction, ray)) != glm::vec3(0, 0, 0)){
                if(spec){
                    output += attenuation * L;
                } else if(prevPDF > 0){
                    float light_pdf = lightSampler->PMF(interaction.AreaLight) * interaction.AreaLight->PDF(interaction, ray);
                    float w = prevPDF * prevPDF / (prevPDF * prevPDF + light_pdf * light_pdf);
                    output += attenuation * L * w;
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


            std::optional<BxDFSample> bxdf = interaction.mat->scatter(ray, interaction, new_ray, scatter_random_variable, random_variables);
            if(!bxdf){
                return output;//absorbed   
            }


            new_ray.medium = interaction.getMedium(new_ray.dir);
            new_ray.time = ray.time;

            if(EnableCaustics && interaction.mat->is_specular(interaction)){
                //&& !bxdf->isSpecular() 
                //is dielectric
                attenuation *= bxdf->f * std::abs(glm::dot(interaction.ns, new_ray.dir)) / bxdf->pdf;

                prevPDF = interaction.mat->PDF(ray, interaction, new_ray);
                spec = false;
                glm::vec3 acc = attenuation * SampleLd(ray, interaction, light_selection_random_variable, light_random_variables);
                output += acc;
                if(acc == glm::vec3(0, 0, 0))spec = true;
                //maybe do same for general reflection not just dielectrics


                ray = new_ray;
                continue;
            }


            //this helps when medium intersect another object
            //when we bounce we will be in same medium as before
            //this breaks if we go into specular surface so do this if also not specular?
            if(!bxdf->isTransmissive() && glm::dot(ray.dir, interaction.ns) <= 0) //transmissive flag
                new_ray.medium = ray.medium;


            spec = bxdf->isSpecular();
            if(!spec){
                output += attenuation * SampleLd(ray, interaction, light_selection_random_variable, light_random_variables);
                prevPDF = interaction.mat->PDF(ray, interaction, new_ray);
            }


            attenuation *= bxdf->f * std::abs(glm::dot(interaction.ns, new_ray.dir)) / bxdf->pdf;

            ray = new_ray;
        }

        if(rr_depth++ > 3){
            float rr_prob = std::fmin(0.95f, std::fmaxf(std::fmaxf(attenuation.x, attenuation.y), attenuation.z));

            if(rr_random_variable >= rr_prob)break;
            attenuation /= rr_prob;
        }

    }

    return output;
}


glm::vec3 VolPathIntegrator::SampleLd(const Ray& ray, const GeometricInteraction& interaction, float u, const glm::vec2& UV) const{
    std::shared_ptr<Light> sampled_light = lightSampler->Sample(u);
    if(sampled_light == nullptr)return { 0,0,0 };
    glm::vec3 Tr = { 1,1,1 };
    SurfaceInteraction intr;

    LightSample lightSample = sampled_light->sample(UV, ray.time);
    glm::vec3 lightDir;

    float t = 0;

    if(lightSample.isDeltaInteraction()){
        lightDir = lightSample.dir;
        t = std::numeric_limits<float>::infinity();
    } else{
        lightDir = lightSample.interaction.p - interaction.p;
        t = glm::length(lightDir) - shadowEpsilon;//was 0.0001f
        t -= shadowEpsilon;
    }

    Ray shadow_ray(interaction.p, glm::normalize(lightDir), ray.time, ray.medium);

    float light_pdf = lightSampler->PMF(sampled_light);
    if(light_pdf <= 0) return { 0,0,0 };

    glm::vec3 f;
    float samplingPDF;
    if(interaction.isMediumInteraction()){
        samplingPDF = static_cast<const MediumInteraction*>(&interaction)->phaseFunction->PDF(ray.dir, shadow_ray.dir);
        f = glm::vec3(samplingPDF);
    } else{
        const SurfaceInteraction* surfIntr = static_cast<const SurfaceInteraction*>(&interaction);
        float dot = glm::dot(surfIntr->ns, shadow_ray.dir);
        if(dot * glm::dot(ray.dir, surfIntr->ns) >= 0)return { 0,0,0 };
        samplingPDF = surfIntr->mat->PDF(ray, *surfIntr, shadow_ray);
        f = surfIntr->mat->calc_attenuation(ray, *surfIntr, shadow_ray) * std::abs(dot);//maybe need abs?
    }


    if(f == glm::vec3(0, 0, 0) || IntersectTr(shadow_ray, intr, Tr, t))return { 0,0,0 };

    if(sampled_light->isDelta()){
        return Tr * lightSample.L * f / light_pdf;
    } else{//if(!lightSample.isDeltaInteraction()) prevents n being (0,0,0) -> good for dengenerate triagnle, bad for infinite area lights
        light_pdf *= sampled_light->PDF(lightSample.interaction, shadow_ray);
        if(light_pdf <= 0)return { 0,0,0 };
        float w2 = light_pdf * light_pdf;
        float w1 = samplingPDF;
        w1 = w1 * w1;
        float w_light = (w2) / (w1 + w2);
        return Tr * sampled_light->L(lightSample.interaction, shadow_ray) * f * w_light / light_pdf;
    }
    return { 0,0,0 };
}

