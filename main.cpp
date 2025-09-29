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
#include "Interaction.hpp"

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
//#include <OpenImageDenoise/oidn.hpp>
#include "Scene.hpp"
#include "ResourceManager.hpp"

#include "Integrators.hpp"

/*
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

    
    double halfWidth  = std::tan(fov * 0.5f);
    double halfHeight = halfWidth * height / width;
    glm::dvec3 up = {0,1,0};
    glm::dvec3 w = glm::normalize(lookfrom-lookat);
    glm::dvec3 u = glm::normalize(glm::cross(up,w));
    glm::dvec3 v = glm::cross(w,u);

    double defocus_radius = focus_dist * std::tan(defocus_angle/2.0f);
    glm::dvec3 defocus_disk_u = u * defocus_radius;
    glm::dvec3 defocus_disk_v = v * defocus_radius;
    
    unsigned int threads = std::thread::hardware_concurrency();
    std::vector<std::thread> workers;
    std::atomic<int> done{0};
 



    int samples = 64;//64*16*4 -> 4 hours
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
                            if(std::isnan(color.x) ||  std::isnan(color.y) ||  std::isnan(color.z)){
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
  
    film->WritePPM(outputImage);
}
*/

void NoModel(){
    auto scene = std::make_shared<Scene>();
    auto white = ResourceManager::get_instance().GetTexture<SolidColor>("whiteTexture",glm::vec3(.9));
    auto green = ResourceManager::get_instance().GetTexture<SolidColor>("greenTexture",glm::vec3(.2,.3,.1));
    auto light = std::make_shared<lambertian>(glm::vec3(0));
    auto ch =  std::make_shared<lambertian>(glm::vec3{.2,.3,.1});
    auto glass = std::make_shared<dielectric>(1.5,glm::vec3(1));
    auto checker = std::make_shared<lambertian>(ResourceManager::get_instance().GetTexture<CheckerTexture>("greenTexture",white,green,glm::vec2{0.02}));
    ch =  std::make_shared<lambertian>(std::make_shared<CheckerTexture>(white,green,glm::vec2{0.001,0.001}));
    std::shared_ptr<AreaLight> area = std::make_shared<AreaLight>(std::make_shared<QuadShape>(glm::vec3(0.3,1.5,0), glm::vec3(-0.15,0,0), glm::vec3(0,0,-0.15)),glm::vec3(600),false);
    auto outsideMedium = std::make_shared<HomogeneusMedium>(glm::vec3{0.01f, 0.9f, 0.9f},glm::vec3{1.0f, 0.1f, 0.1f},std::make_shared<HenyeyGreenstein>(0.9),0.5f);


    scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<QuadShape>(glm::vec3(0.3,1.5,0), glm::vec3(-0.15,0,0), glm::vec3(0,0,-0.15)), light, area, nullptr));//-0.3, -1
    scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<QuadShape>(glm::vec3(-100,-0.3,-100), glm::vec3(1000,0,0), glm::vec3(0,0,1000)), ch, nullptr, nullptr));
    //scene->Add(new Model("/home/markov/Documents/Coding/CPP/raytracing_in_one_weekend/temp_other.assbin"));


    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(0,0.1,-1.2),0.5),std::make_shared<lambertian>(glm::vec3(0.1, 0.2, 0.5)),nullptr));
    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(-1,0,-1),0.5),std::make_shared<dielectric>(1.5),nullptr));
    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(-1,0,-1),0.4),std::make_shared<dielectric>(1/1.5),nullptr));
    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(1,0,-1),0.5),std::make_shared<metal>(glm::vec3(0.8, 0.6, 0.2)),nullptr));
    
    scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<SphereShape>(glm::vec3(1,0,-1),0.5),nullptr,nullptr,std::make_shared<HomogeneusMedium>(glm::vec3{0.01f, 0.9f, 0.9f},glm::vec3{1.0f, 0.1f, 0.1f},std::make_shared<HenyeyGreenstein>(0.8),25.0f,glm::vec3{1,1,1},0)));//Was 1 
        //d_list[1] = new Sphere{glm::vec3(-0.8,1,-0.5), 0.5,
        //                        new Light(glm::vec3(8, 8, 8))};
    //scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<SphereShape>(glm::vec3(0,0,0),20),nullptr,nullptr,outsideMedium));
    

    std::shared_ptr<LightSampler> ls = std::make_shared<PowerLightSampler>();
    scene->PreProcess();
    ls->Add(scene->GetLights());
    //ls->Add(std::make_shared<PointLight>(glm::vec3(0.3,1.5,0),glm::vec3(6)));
    ls->PreProcess(scene->BoundingBox());
    
    double fov = 1.7;


    glm::dvec3 lookfrom = {-1000,300,0};
    lookfrom = {17.3,1.2,7.2};

    lookfrom = {0.3,0.4,1};//dragon

    glm::dvec3 lookat = {0,300,0};
    lookat = {0,0,0};

   

    int samples = 16;



    int sqrts = std::sqrt(samples);

    std::shared_ptr<Film> film = std::make_shared<Film>(glm::ivec2{1920,1080},std::make_shared<MitchellFilter>());

    
    auto camera = std::make_shared<Camera>(lookfrom,lookat,fov,film);
    auto sampler = std::make_shared<StratifiedSampler>(sqrts,sqrts);
    auto integrator = std::make_shared<VolPathIntegrator>(scene,camera,sampler,ls,128);

    //camera->SetMedium(outsideMedium);
    //scene->SetMedium(outsideMedium);

    integrator->Render();
    camera->GetFilm()->WritePNG("RenderedScene");
    camera->GetFilm()->WritePPM("RenderedScene");
    camera->GetFilm()->WriteJPG("RenderedScene",100);
    ResourceManager::get_instance().releaseTextures();
}

void Miguel(){
    auto scene = std::make_shared<Scene>();
    auto white = ResourceManager::get_instance().GetTexture<SolidColor>("whiteTexture",glm::vec3(.9));
    auto green = ResourceManager::get_instance().GetTexture<SolidColor>("greenTexture",glm::vec3(.2,.3,.1));
    auto light = std::make_shared<lambertian>(glm::vec3(0));
    auto ch =  std::make_shared<lambertian>(glm::vec3{.2,.3,.1});
    auto glass = std::make_shared<dielectric>(1.5,glm::vec3(1));
    auto checker = std::make_shared<lambertian>(ResourceManager::get_instance().GetTexture<CheckerTexture>("greenTexture",white,green,glm::vec2{0.02}));
    ch =  std::make_shared<lambertian>(std::make_shared<CheckerTexture>(white,green,glm::vec2{0.001,0.001}));
    std::shared_ptr<AreaLight> area = std::make_shared<AreaLight>(std::make_shared<QuadShape>(glm::vec3(0.3,1.5,0), glm::vec3(-0.15,0,0), glm::vec3(0,0,-0.15)),glm::vec3(600),false);

    scene->Add(ResourceManager::get_instance().GetModel("San Miguel","/home/markov/Documents/Coding/CPP/testing/models/HARD/temp.assbin"));


    std::shared_ptr<LightSampler> ls = std::make_shared<PowerLightSampler>();
    scene->PreProcess();
    ls->Add(scene->GetLights());
    ls->Add(std::make_shared<InfiniteLight>(glm::vec3(-1,6,1),25.f*glm::vec3(1,0.93,0.83)));
    //ls->Add(std::make_shared<PointLight>(glm::vec3(0.3,1.5,0),glm::vec3(6)));
    ls->PreProcess(scene->BoundingBox());
    
    double fov = 1.7;


    glm::dvec3 lookfrom = {-1000,300,0};
    lookfrom = {17.3,1.2,7.2};

    glm::dvec3 lookat = {0,300,0};
    lookat = {0,0,0};

    //64*4
    int samples = 64;//64*16*4 -> 4 hours
    int sqrts = std::sqrt(samples);

    std::shared_ptr<Film> film = std::make_shared<Film>(glm::ivec2{1920,1080},std::make_shared<MitchellFilter>());

    
    auto camera = std::make_shared<Camera>(lookfrom,lookat,fov,film);
    auto sampler = std::make_shared<StratifiedSampler>(sqrts,sqrts);
    auto integrator = std::make_shared<PathIntegrator>(scene,camera,sampler,ls,128);

    //camera->SetMedium(outsideMedium);
    //scene->SetMedium(outsideMedium);

    //renderPrim2(scene,1920,1080,ls,"RenderedScene.ppm");
    integrator->Render();
    camera->GetFilm()->WritePNG("RenderedScene");
    ResourceManager::get_instance().releaseTextures();
}

void temp(){
    auto scene = std::make_shared<Scene>();
    auto white = ResourceManager::get_instance().GetTexture<SolidColor>("whiteTexture",glm::vec3(.9));
    auto green = ResourceManager::get_instance().GetTexture<SolidColor>("greenTexture",glm::vec3(.2,.3,.1));
    auto light = std::make_shared<lambertian>(glm::vec3(0));
    auto ch =  std::make_shared<lambertian>(glm::vec3{.2,.3,.1});
    auto glass = std::make_shared<dielectric>(1.5,glm::vec3(1));
    auto checker = std::make_shared<lambertian>(ResourceManager::get_instance().GetTexture<CheckerTexture>("greenTexture",white,green,glm::vec2{0.02}));
    ch =  std::make_shared<lambertian>(std::make_shared<CheckerTexture>(white,green,glm::vec2{0.001,0.001}));
    std::shared_ptr<AreaLight> area = std::make_shared<AreaLight>(std::make_shared<QuadShape>(glm::vec3(0.3,1.5,0), glm::vec3(-0.15,0,0), glm::vec3(0,0,-0.15)),glm::vec3(600),false);
    auto outsideMedium = std::make_shared<HomogeneusMedium>(glm::vec3{0.01f, 0.9f, 0.9f},glm::vec3{1.0f, 0.1f, 0.1f},std::make_shared<HenyeyGreenstein>(0.9),0.5f);


    scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<QuadShape>(glm::vec3(0.3,1.5,0), glm::vec3(-0.15,0,0), glm::vec3(0,0,-0.15)), light, area, nullptr));//-0.3, -1
    scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<QuadShape>(glm::vec3(-100,-0.3,-100), glm::vec3(1000,0,0), glm::vec3(0,0,1000)), ch, nullptr, nullptr));
    //scene->Add(new Model("/home/markov/Documents/Coding/CPP/raytracing_in_one_weekend/temp_other.assbin"));
    glm::mat4 pos = glm::mat4(1);
    pos = glm::translate(pos,{-0.1,0,-0.1});
    pos = glm::rotate(pos,glm::radians(-60.0f),glm::normalize(glm::vec3(0,1,1)));
    pos = glm::scale(pos,glm::vec3(1,1,1));
    std::shared_ptr<Primitive> transformedModel = std::make_shared<TransformedPrimitive>(ResourceManager::get_instance().GetModel("Glass Dragon","/home/markov/Documents/Coding/CPP/testing/models/temp_other.assbin",glass,std::make_shared<HomogeneusMedium>(glm::vec3{0.01f, 0.9f, 0.9f},glm::vec3{1.0f, 0.1f, 0.1f},std::make_shared<HenyeyGreenstein>(0.8),25.0f)),pos);
  
    //scene->Add(transformedModel);

    scene->Add(ResourceManager::get_instance().GetModel("Medium Dragon","/home/markov/Documents/Coding/CPP/testing/models/temp_other.assbin",
            nullptr,
            std::make_shared<HomogeneusMedium>( glm::vec3{0.01f, 0.9f, 0.9f},
                                                glm::vec3{1.0f, 0.1f, 0.1f},
                                                std::make_shared<HenyeyGreenstein>(0.8),
                                                25.0f)));//was1.2

    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(0,0.1,-1.2),0.5),std::make_shared<lambertian>(glm::vec3(0.1, 0.2, 0.5)),nullptr));
    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(-1,0,-1),0.5),std::make_shared<dielectric>(1.5),nullptr));
    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(-1,0,-1),0.4),std::make_shared<dielectric>(1/1.5),nullptr));
    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(1,0,-1),0.5),std::make_shared<metal>(glm::vec3(0.8, 0.6, 0.2)),nullptr));
    
    scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<SphereShape>(glm::vec3(1,0,-1),0.5),nullptr,nullptr,std::make_shared<HomogeneusMedium>(glm::vec3{0.01f, 0.9f, 0.9f},glm::vec3{1.0f, 0.1f, 0.1f},std::make_shared<HenyeyGreenstein>(0.8),25.0f)));//Was 1 
        //d_list[1] = new Sphere{glm::vec3(-0.8,1,-0.5), 0.5,
        //                        new Light(glm::vec3(8, 8, 8))};
    //scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<SphereShape>(glm::vec3(0,0,0),20),nullptr,nullptr,outsideMedium));
    

    std::shared_ptr<LightSampler> ls = std::make_shared<PowerLightSampler>();
    scene->PreProcess();
    ls->Add(scene->GetLights());
    //ls->Add(std::make_shared<PointLight>(glm::vec3(0.3,1.5,0),glm::vec3(6)));
    ls->PreProcess(scene->BoundingBox());
    
    double fov = 1.7;


    glm::dvec3 lookfrom = {-1000,300,0};
    lookfrom = {17.3,1.2,7.2};

    lookfrom = {0.3,0.4,1};//dragon

    glm::dvec3 lookat = {0,300,0};
    lookat = {0,0,0};

    //64*4
    int samples = 16;//64*16*4 -> 4 hours
    int sqrts = std::sqrt(samples);

    std::shared_ptr<Film> film = std::make_shared<Film>(glm::ivec2{1920,1080},std::make_shared<MitchellFilter>());

    
    auto camera = std::make_shared<Camera>(lookfrom,lookat,fov,film);
    auto sampler = std::make_shared<StratifiedSampler>(sqrts,sqrts);
    auto integrator = std::make_shared<VolPathIntegrator>(scene,camera,sampler,ls,128);

    //camera->SetMedium(outsideMedium);
    //scene->SetMedium(outsideMedium);

    integrator->Render();
    camera->GetFilm()->WritePNG("RenderedScene");
    camera->GetFilm()->WritePPM("RenderedScene");
    camera->GetFilm()->WriteJPG("RenderedScene",100);
    ResourceManager::get_instance().releaseTextures();
}

//how to deal when we spawn ray inside objects with medium?
//interaction spawnRay() -> this spawns ray but gives us 
// ### />   <- .
//getMedium would be ###, but we would set medium to be null (coming from outside)
// </###   <- .
//getMedium would be nullptr, but we would set medium to be ### (coming from inside)

int main(){
    stbi_set_flip_vertically_on_load(true);
    //"/home/markov/Documents/Coding/CPP/testing/stanford/common-3d-test-models-master/data/lucy.obj"
    switch(0){
        case 0:
            //wierd new fireflys?
            temp();
            break;
        case 1:
            Miguel();
            break;
        case 2:
            NoModel();
        break;
    }
    ResourceManager::get_instance().releaseTextures();
}
