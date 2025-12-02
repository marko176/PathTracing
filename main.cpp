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



void MatTest(){
    auto scene = std::make_shared<Scene>();
    auto white = ResourceManager::get_instance().GetTexture<SolidColor>("whiteTexture", glm::vec3(.9));
    auto green = ResourceManager::get_instance().GetTexture<SolidColor>("greenTexture", glm::vec3(.2, .3, .1));
    auto light = std::make_shared<MicrofacetDiffuse>(glm::vec3(0));
    auto ch = std::make_shared<MicrofacetDiffuse>(glm::vec3 { .2,.3,.1 });
    auto glass = std::make_shared<MicrofacetDielectric>(1.5, 0.0, glm::vec3(1));
    auto checker = std::make_shared<MicrofacetDiffuse>(ResourceManager::get_instance().GetTexture<CheckerTexture>("greenTexture", white, green, glm::vec2 { 0.02 }));
    ch = std::make_shared<MicrofacetDiffuse>(std::make_shared<CheckerTexture>(white, green, glm::vec2 { 0.001,0.001 }));
    std::shared_ptr<AreaLight> area = std::make_shared<AreaLight>(std::make_shared<QuadShape>(glm::vec3(0.3, 1.5, 0), glm::vec3(-0.15, 0, 0), glm::vec3(0, 0, -0.15)), glm::vec3(600), false);
    auto outsideMedium = std::make_shared<HomogeneusMedium>(glm::vec3 { 0.01f, 0.9f, 0.9f }, glm::vec3 { 1.0f, 0.1f, 0.1f }, std::make_shared<HenyeyGreenstein>(0.9), 0.5f);



    scene->Add(std::make_shared<GeometricPrimitive>(area->getShape(), light, area, nullptr));//-0.3, -1
    scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<QuadShape>(glm::vec3(-100, -0.3, -100), glm::vec3(1000, 0, 0), glm::vec3(0, 0, 1000)), ch, nullptr, nullptr));
    auto b = std::make_shared<SolidColor>(glm::vec3 { .1,.2,.5 });
    std::vector<float> values = { 1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.0005f,0.0,0.05,0.01,0.005,0.001,0.0001,0.000 };
    std::vector<float> metalValues = { 0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1 };
    for(std::size_t i = 0;i < values.size();i++){
        for(std::size_t j = 0;j < metalValues.size();j++){
            auto mat = std::make_shared<MicrofacetDiffuse>(b, nullptr, std::make_shared<SolidColor>(glm::vec3(values[i])), std::make_shared<SolidColor>(glm::vec3(metalValues[j])));
            auto matee = std::make_shared<SpecularConductor>(glm::vec3 { .1,.2,.5 });

            scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<SphereShape>(glm::vec3(i - (int)values.size() / 2, j + 0.5, 0), 0.4), mat));
        }
    }
    auto mat1 = std::make_shared<MicrofacetDiffuse>(b, nullptr, std::make_shared<SolidColor>(glm::vec3(1)), std::make_shared<SolidColor>(glm::vec3(0)));
    auto mat2 = std::make_shared<MicrofacetDiffuse>(b, nullptr, std::make_shared<SolidColor>(glm::vec3(0.7)), std::make_shared<SolidColor>(glm::vec3(0)));
    auto mat3 = std::make_shared<MicrofacetDiffuse>(b, nullptr, std::make_shared<SolidColor>(glm::vec3(0.5)), std::make_shared<SolidColor>(glm::vec3(0)));
    auto mat4 = std::make_shared<MicrofacetDiffuse>(b, nullptr, std::make_shared<SolidColor>(glm::vec3(0.3)), std::make_shared<SolidColor>(glm::vec3(0)));
    auto mat5 = std::make_shared<MicrofacetDiffuse>(b, nullptr, std::make_shared<SolidColor>(glm::vec3(0.1)), std::make_shared<SolidColor>(glm::vec3(0)));
    auto mat6 = std::make_shared<MicrofacetDiffuse>(b, nullptr, std::make_shared<SolidColor>(glm::vec3(0.05)), std::make_shared<SolidColor>(glm::vec3(0)));
    auto mat7 = std::make_shared<MicrofacetDiffuse>(b, nullptr, std::make_shared<SolidColor>(glm::vec3(0.0011)), std::make_shared<SolidColor>(glm::vec3(0)));
    auto mat8 = std::make_shared<MicrofacetDiffuse>(b, nullptr, std::make_shared<SolidColor>(glm::vec3(0.0006)), std::make_shared<SolidColor>(glm::vec3(0)));
    auto mat9 = std::make_shared<MicrofacetDiffuse>(b, nullptr, std::make_shared<SolidColor>(glm::vec3(0.0001)), std::make_shared<SolidColor>(glm::vec3(0)));

    auto matall = std::make_shared<SolidColor>(glm::vec3(1));
    auto met1 = std::make_shared<MicrofacetDiffuse>(b, nullptr, std::make_shared<SolidColor>(glm::vec3(1)), matall);
    auto met2 = std::make_shared<MicrofacetDiffuse>(b, nullptr, std::make_shared<SolidColor>(glm::vec3(0.7)), matall);
    auto met3 = std::make_shared<MicrofacetDiffuse>(b, nullptr, std::make_shared<SolidColor>(glm::vec3(0.5)), matall);
    auto met4 = std::make_shared<MicrofacetDiffuse>(b, nullptr, std::make_shared<SolidColor>(glm::vec3(0.3)), matall);
    auto met5 = std::make_shared<MicrofacetDiffuse>(b, nullptr, std::make_shared<SolidColor>(glm::vec3(0.1)), matall);
    auto met6 = std::make_shared<MicrofacetDiffuse>(b, nullptr, std::make_shared<SolidColor>(glm::vec3(0.05)), matall);
    auto met7 = std::make_shared<MicrofacetDiffuse>(b, nullptr, std::make_shared<SolidColor>(glm::vec3(0.01)), matall);
    auto met8 = std::make_shared<MicrofacetDiffuse>(b, nullptr, std::make_shared<SolidColor>(glm::vec3(0.005)), matall);
    auto met9 = std::make_shared<MicrofacetDiffuse>(b, nullptr, std::make_shared<SolidColor>(glm::vec3(0.0001)), matall);
    std::shared_ptr<AreaLight> animatedArea = std::make_shared<AreaLight>(std::make_shared<SphereShape>(glm::vec3(0, 0.1, -1.2), 0.5), glm::vec3(10), false);
    auto tm = std::make_shared<GeometricPrimitive>(std::make_shared<SphereShape>(glm::vec3(0, 0.1, -1.2), 0.5), mat9);
    //auto animeted = std::make_shared<AnimetedPrimitive>(tm,glm::vec3{0,0.2,0},glm::vec2{0,1});
    //scene->Add(tm);

    //scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<SphereShape>(glm::vec3(-1,0,-1),0.5),mat8));
    //scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<SphereShape>(glm::vec3(-1,0,-1),0.4),std::make_shared<dielectric>(1/1.5),nullptr));
    //scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<SphereShape>(glm::vec3(-1,0,0.2),0.5),mat5));

    //scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<SphereShape>(glm::vec3(1,0,-1),0.5),mat4));//Was 1 
        //d_list[1] = new Sphere{glm::vec3(-0.8,1,-0.5), 0.5,
        //                        new Light(glm::vec3(8, 8, 8))};
    //scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<SphereShape>(glm::vec3(0,0,0),20),nullptr,nullptr,outsideMedium));

    auto lightFunc = [](const Ray& ray){
        float a = 0.5f * (ray.dir.y + 1.0f);
        return 1.5f * ((1.0f - a) * glm::vec3(1, 0.85, 0.55) + a * glm::vec3(0.45, 0.65, 1));
        };

    //scene->infiniteLights.push_back(std::make_shared<UniformInfiniteLight>(glm::vec3{0,0,1}));
    scene->infiniteLights.push_back(std::make_shared<FunctionInfiniteLight>(lightFunc));
    //scene->infiniteLights.push_back(std::make_shared<TextureInfiniteLight>(std::make_shared<FloatImageTexture>("/home/markov/Downloads/kloofendal_48d_partly_cloudy_puresky_8k.hdr"),600,[](float r){return 4 * std::sqrt(r);}));
    //scene->infiniteLights.push_back(std::make_shared<TextureInfiniteLight>(std::make_shared<FloatImageTexture>("/home/markov/Downloads/lilienstein_4k.hdr"),600,[](float r){return 4 * std::sqrt(r);}));

    std::shared_ptr<LightSampler> ls = std::make_shared<PowerLightSampler>();
    scene->BuildTlas<TLAS4>();


    ls->Add(scene->GetLights());
    //ls->Add(std::make_shared<PointLight>(glm::vec3(0.3,1.5,0),glm::vec3(6)));
    ls->PreProcess(scene->BoundingBox());

    double fov = 1.65;


    glm::dvec3 lookfrom = { -1000,300,0 };
    lookfrom = { 17.3,1.2,7.2 };

    lookfrom = { 0,metalValues.size() / 2 + 0.5,9 };//dragon

    glm::dvec3 lookat = { 0,300,0 };
    lookat = { 0,metalValues.size() / 2 + 0.5,0 };



    int samples = 100;



    int sqrts = std::sqrt(samples);

    std::shared_ptr<Film> film = std::make_shared<Film>(glm::ivec2 { 1920,1080 }, std::make_shared<MitchellFilter>());


    auto camera = std::make_shared<Camera>(lookfrom, lookat, fov, film, glm::vec2(0, 1));
    auto sampler = std::make_shared<StratifiedSampler>(sqrts, sqrts);
    auto integrator = std::make_shared<VolPathIntegrator>(scene, camera, sampler, ls, 128);

    //camera->SetMedium(outsideMedium);
    //scene->SetMedium(outsideMedium);

    integrator->Render();
    camera->GetFilm()->WritePNG("RenderedScene");
    camera->GetFilm()->WritePPM("RenderedScene");
    camera->GetFilm()->WriteJPG("RenderedScene", 90);
    ResourceManager::get_instance().releaseTextureCache();
    ResourceManager::get_instance().releaseMeshCache();
    ResourceManager::get_instance().releaseModelCache();
}

void NoModel(){
    auto scene = std::make_shared<Scene>();
    auto white = ResourceManager::get_instance().GetTexture<SolidColor>("whiteTexture", glm::vec3(.9));
    auto green = ResourceManager::get_instance().GetTexture<SolidColor>("greenTexture", glm::vec3(.2, .3, .1));
    auto light = std::make_shared<MicrofacetDiffuse>(glm::vec3(0));
    auto ch = std::make_shared<MicrofacetDiffuse>(glm::vec3 { .2,.3,.1 });
    auto glass = std::make_shared<MicrofacetDielectric>(1.5, 0.0, glm::vec3(1));
    auto checker = std::make_shared<MicrofacetDiffuse>(ResourceManager::get_instance().GetTexture<CheckerTexture>("greenTexture", white, green, glm::vec2 { 0.02 }));
    ch = std::make_shared<MicrofacetDiffuse>(std::make_shared<CheckerTexture>(white, green, glm::vec2 { 0.001,0.001 }));
    std::shared_ptr<AreaLight> area = std::make_shared<AreaLight>(std::make_shared<QuadShape>(glm::vec3(0.3, 2.5, 0), glm::vec3(-0.15, 0, 0), glm::vec3(0, 0, -0.15)), glm::vec3(1000), false);
    auto outsideMedium = std::make_shared<HomogeneusMedium>(glm::vec3 { 0.01f, 0.9f, 0.9f }, glm::vec3 { 1.0f, 0.1f, 0.1f }, std::make_shared<HenyeyGreenstein>(0.75), 0.1f);


    scene->Add(std::make_shared<GeometricPrimitive>(area->getShape(), light, area, nullptr));//-0.3, -1


    auto micro = std::make_shared<MicrofacetDielectric>(1.5, 0.15, glm::vec3 { 1 });
    scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<QuadShape>(glm::vec3(-100, -0.3, -100), glm::vec3(1000, 0, 0), glm::vec3(0, 0, 1000)), ch, nullptr, nullptr));




    //scene->Add(new Model("/home/markov/Documents/Coding/CPP/raytracing_in_one_weekend/temp_other.assbin"));
    std::shared_ptr<AreaLight> animatedArea = std::make_shared<AreaLight>(std::make_shared<SphereShape>(glm::vec3(0, 0.1, -1.2), 0.5), glm::vec3(10), false);
    auto tm = std::make_shared<GeometricPrimitive>(std::make_shared<SphereShape>(glm::vec3(0, 0.1, -1.2), 0.5), std::make_shared<MicrofacetDiffuse>(glm::vec3(0.1, 0.2, 0.5)), nullptr);
    auto animated = std::make_shared<AnimatedPrimitive>(tm, glm::vec3 { 0,0.2,0 }, glm::vec2 { 0,1 });
    scene->Add(animated);

    //scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<SphereShape>(glm::vec3(-1,0,-1),0.5),std::make_shared<dielectric>(1.5),nullptr));



    scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<SphereShape>(glm::vec3(-1, 0.3, -1), 0.5), micro, nullptr));





    //scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<SphereShape>(glm::vec3(-1,0,-1),0.4),std::make_shared<dielectric>(1/1.5),nullptr));

    auto met = std::make_shared<MicrofacetDiffuse>(std::make_shared<SolidColor>(glm::vec3(0.8, 0.6, 0.2)), nullptr, std::make_shared<SolidColor>(glm::vec3(0)), std::make_shared<SolidColor>(glm::vec3(1)));
    auto metRef = std::make_shared<SpecularConductor>(glm::vec3(0.8, 0.6, 0.2));
    scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<SphereShape>(glm::vec3(-1, 0, 0.2), 0.5), met, nullptr));

    scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<SphereShape>(glm::vec3(1, 0, -1), 0.5), nullptr, nullptr, std::make_shared<HomogeneusMedium>(glm::vec3 { 0.01f, 0.9f, 0.9f }, glm::vec3 { 1.0f, 0.1f, 0.1f }, std::make_shared<HenyeyGreenstein>(0.8), 5.0f, glm::vec3 { 1,1,1 }, 0)));//Was 1 
    //d_list[1] = new Sphere{glm::vec3(-0.8,1,-0.5), 0.5,
    //                        new Light(glm::vec3(8, 8, 8))};
//scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<SphereShape>(glm::vec3(0,0,0),20),nullptr,nullptr,outsideMedium));

    
    auto lightFunc = [](const Ray& ray){
        float a = 0.5f * (ray.dir.y + 1.0f);
        return 3.f * ((1.0f - a) * glm::vec3(1, 0.85, 0.55) + a * glm::vec3(0.45, 0.65, 1));
        };
    

    //scene->infiniteLights.push_back(std::make_shared<UniformInfiniteLight>(glm::vec3{0,0,1}));
    //scene->infiniteLights.push_back(std::make_shared<FunctionInfiniteLight>(lightFunc));
    //scene->infiniteLights.push_back(std::make_shared<TextureInfiniteLight>(std::make_shared<FloatImageTexture>("/home/markov/Downloads/kloofendal_48d_partly_cloudy_puresky_8k.hdr"), 600 / 255.0f, [](float r){return 4 * std::sqrt(r);}));
    //scene->infiniteLights.push_back(std::make_shared<TextureInfiniteLight>(std::make_shared<FloatImageTexture>("/home/markov/Downloads/lilienstein_4k.hdr"),600/255.0f,[](float r){return 4 * std::sqrt(r);}));
    //scene->infiniteLights.push_back(std::make_shared<TextureInfiniteLight>(std::make_shared<FloatImageTexture>("/home/markov/Downloads/shanghai_bund_8k.hdr"),600/255.0f,[](float r){return 4 * std::sqrt(r);}));



    std::shared_ptr<LightSampler> ls = std::make_shared<PowerLightSampler>();

    scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<SphereShape>(glm::vec3(0, 0, 0), 10), checker, nullptr));


    scene->BuildTlas<TLAS4>();


    ls->Add(scene->GetLights());
    //ls->Add(std::make_shared<PointLight>(glm::vec3(0.3,1.5,0),glm::vec3(6)));
    ls->PreProcess(scene->BoundingBox());

    double fov = 2.5;//1.7 , 2.5


    glm::dvec3 lookfrom = { -1000,300,0 };
    lookfrom = { 17.3,1.2,7.2 };

    lookfrom = { 0.3,0.4,1 };//dragon
    lookfrom = { 0.3,1.7,1 };
    glm::dvec3 lookat = { 0,300,0 };
    lookat = { 0,1.5,0 };



    int samples = 100*4;



    int sqrts = std::sqrt(samples);

    std::shared_ptr<Film> film = std::make_shared<Film>(glm::ivec2 { 1920,1080 }, std::make_shared<MitchellFilter>());


    auto camera = std::make_shared<Camera>(lookfrom, lookat, fov, film, glm::vec2(0, 1));
    auto sampler = std::make_shared<StratifiedSampler>(sqrts, sqrts);
    auto integrator = std::make_shared<VolPathIntegrator>(scene, camera, sampler, ls, 64);

    camera->SetMedium(outsideMedium);
    //scene->SetMediumBounds({-5,-5,-5},{5,5,5});this should be primitive which we can also hit
    scene->SetMedium(outsideMedium);

    integrator->Render();
    camera->GetFilm()->WritePNG("RenderedScene");
    camera->GetFilm()->WritePPM("RenderedScene");
    camera->GetFilm()->WriteJPG("RenderedScene", 90);
    ResourceManager::get_instance().releaseTextureCache();
    ResourceManager::get_instance().releaseMeshCache();
    ResourceManager::get_instance().releaseModelCache();
}

void Miguel(){
    auto scene = std::make_shared<Scene>();
    auto white = ResourceManager::get_instance().GetTexture<SolidColor>("whiteTexture", glm::vec3(.9));
    auto green = ResourceManager::get_instance().GetTexture<SolidColor>("greenTexture", glm::vec3(.2, .3, .1));
    auto light = std::make_shared<MicrofacetDiffuse>(glm::vec3(0));
    auto ch = std::make_shared<MicrofacetDiffuse>(glm::vec3 { .2,.3,.1 });
    auto glass = std::make_shared<MicrofacetDielectric>(1.5, 0.0, glm::vec3(1));
    auto checker = std::make_shared<MicrofacetDiffuse>(ResourceManager::get_instance().GetTexture<CheckerTexture>("greenTexture", white, green, glm::vec2 { 0.02 }));
    ch = std::make_shared<MicrofacetDiffuse>(std::make_shared<CheckerTexture>(white, green, glm::vec2 { 0.001,0.001 }));
    std::shared_ptr<AreaLight> area = std::make_shared<AreaLight>(std::make_shared<QuadShape>(glm::vec3(0.3, 1.5, 0), glm::vec3(-0.15, 0, 0), glm::vec3(0, 0, -0.15)), glm::vec3(600), false);

    scene->Add(ResourceManager::get_instance().CacheModel<BLAS4>("San Miguel", "/home/markov/Documents/Coding/CPP/testing/models/HARD/temp.assbin"));

    auto lightFunc = [](const Ray& ray){
        float a = 0.5f * (ray.dir.y + 1.0f);
        return 1.5f * ((1.0f - a) * glm::vec3(1, 0.85, 0.55) + a * glm::vec3(0.45, 0.65, 1));
        };

    //scene->infiniteLights.push_back(std::make_shared<UniformInfiniteLight>(glm::vec3{0,0,1}));
    scene->infiniteLights.push_back(std::make_shared<FunctionInfiniteLight>(lightFunc));
    //scene->infiniteLights.push_back(std::make_shared<TextureInfiniteLight>(std::make_shared<FloatImageTexture>("/home/markov/Downloads/noon_grass_8k.hdr"),600/255.0f,[](float r){return 4 * std::sqrt(r);}));

    std::shared_ptr<LightSampler> ls = std::make_shared<PowerLightSampler>();
    scene->BuildTlas<TLAS4>();
    ls->Add(scene->GetLights());
    ls->Add(std::make_shared<DistantLight>(glm::vec3(-1, 6, 1), 25.f * glm::vec3(1, 0.93, 0.83)));
    //ls->Add(std::make_shared<PointLight>(glm::vec3(0.3,1.5,0),glm::vec3(6)));
    ls->PreProcess(scene->BoundingBox());

    double fov = 1.7;


    glm::dvec3 lookfrom = { -1000,300,0 };
    lookfrom = { 17.3,1.2,7.2 };

    glm::dvec3 lookat = { 0,300,0 };
    lookat = { 0,0,0 };

    //64*4
    //100 -> ~618 s
    //2k 100 -> 1100 s
    //100 -> 572 s with bxdf sample
    //8% improvement
    //100 -> 516 simd
    //100 -> 550 simd + emissive
    //1024 at 1080p 5027040ms

    //567 bvh2
    //518 simd bvh2
    //384s BVH4 SIMD 386
    //367s BVH8 SIMD 372

    //BVH4 363 -> randomAVXx8
    //BVH8 349 -> randomAVXx8 
    int samples = 100;
    int sqrts = std::sqrt(samples);

    std::shared_ptr<Film> film = std::make_shared<Film>(glm::ivec2 { 1920,1080 }, std::make_shared<MitchellFilter>());


    auto camera = std::make_shared<Camera>(lookfrom, lookat, fov, film);
    auto sampler = std::make_shared<StratifiedSampler>(sqrts, sqrts);
    auto integrator = std::make_shared<PathIntegrator>(scene, camera, sampler, ls, 64);

    //camera->SetMedium(outsideMedium);
    //scene->SetMedium(outsideMedium);

    //renderPrim2(scene,1920,1080,ls,"RenderedScene.ppm");
    integrator->Render();
    camera->GetFilm()->WritePPM("RenderedScene");
    camera->GetFilm()->WritePNG("RenderedScene");
    camera->GetFilm()->WriteJPG("RenderedScene", 100);
    ResourceManager::get_instance().releaseTextureCache();
    ResourceManager::get_instance().releaseMeshCache();
    ResourceManager::get_instance().releaseModelCache();
}

void temp(){
    auto scene = std::make_shared<Scene>();
    auto white = ResourceManager::get_instance().GetTexture<SolidColor>("whiteTexture", glm::vec3(.9));
    auto green = ResourceManager::get_instance().GetTexture<SolidColor>("greenTexture", glm::vec3(.2, .3, .1));
    auto light = std::make_shared<MicrofacetDiffuse>(glm::vec3(0));
    auto ch = std::make_shared<MicrofacetDiffuse>(glm::vec3 { .2,.3,.1 });
    auto glass = std::make_shared<MicrofacetDielectric>(1.5, 0.0, glm::vec3(1));
    auto checker = std::make_shared<MicrofacetDiffuse>(ResourceManager::get_instance().GetTexture<CheckerTexture>("greenTexture", white, green, glm::vec2 { 0.02 }));
    ch = std::make_shared<MicrofacetDiffuse>(std::make_shared<CheckerTexture>(white, green, glm::vec2 { 0.001,0.001 }));
    std::shared_ptr<AreaLight> area = std::make_shared<AreaLight>(std::make_shared<QuadShape>(glm::vec3(0.3, 1.5, 0), glm::vec3(-0.15, 0, 0), glm::vec3(0, 0, -0.15)), glm::vec3(600), false);
    auto outsideMedium = std::make_shared<HomogeneusMedium>(glm::vec3 { 0.01f, 0.9f, 0.9f }, glm::vec3 { 1.0f, 0.1f, 0.1f }, std::make_shared<HenyeyGreenstein>(0.9), 0.5f);


    scene->Add(std::make_shared<GeometricPrimitive>(area->getShape(), light, area, nullptr));//-0.3, -1
    scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<QuadShape>(glm::vec3(-100, -0.3, -100), glm::vec3(1000, 0, 0), glm::vec3(0, 0, 1000)), ch, nullptr, nullptr));
    //scene->Add(new Model("/home/markov/Documents/Coding/CPP/raytracing_in_one_weekend/temp_other.assbin"));
    glm::mat4 pos = glm::mat4(1);
    pos = glm::translate(pos, { -0.1,0,-0.1 });
    pos = glm::rotate(pos, glm::radians(-60.0f), glm::normalize(glm::vec3(0, 1, 1)));
    pos = glm::scale(pos, glm::vec3(1, 1, 1));
    std::shared_ptr<Primitive> transformedModel = std::make_shared<TransformedPrimitive>(ResourceManager::get_instance().CacheModel<BLAS4>("Glass Dragon", "/home/markov/Documents/Coding/CPP/testing/models/temp_other.assbin", glass, std::make_shared<HomogeneusMedium>(glm::vec3 { 0.01f, 0.9f, 0.9f }, glm::vec3 { 1.0f, 0.1f, 0.1f }, std::make_shared<HenyeyGreenstein>(0.8), 25.0f)), pos);

    //scene->Add(transformedModel);
    auto micro = std::make_shared<MicrofacetDielectric>(1.5, 0.0, glm::vec3 { 1 });

    scene->Add(ResourceManager::get_instance().CacheModel<BLAS4>("Medium Dragon", "/home/markov/Documents/Coding/CPP/testing/models/temp_other.assbin",
        ch,
        std::make_shared<HomogeneusMedium>(glm::vec3 { 0.01f, 0.9f, 0.9f },
            glm::vec3 { 1.0f, 0.1f, 0.1f },
            std::make_shared<HenyeyGreenstein>(0.8),
            0.0f,//was 25
            glm::vec3 { 1,1,1 },
            0)));//was1.2

    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(0,0.1,-1.2),0.5),std::make_shared<MicrofacetDiffuse>(glm::vec3(0.1, 0.2, 0.5)),nullptr));
    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(-1,0,-1),0.5),std::make_shared<dielectric>(1.5),nullptr));
    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(-1,0,-1),0.4),std::make_shared<dielectric>(1/1.5),nullptr));
    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(1,0,-1),0.5),std::make_shared<metal>(glm::vec3(0.8, 0.6, 0.2)),nullptr));

    //scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<SphereShape>(glm::vec3(1,0,-1),0.5),nullptr,nullptr,std::make_shared<HomogeneusMedium>(glm::vec3{0.01f, 0.9f, 0.9f},glm::vec3{1.0f, 0.1f, 0.1f},std::make_shared<HenyeyGreenstein>(0.8),25.0f,glm::vec3{1,1,1},1)));//Was 1 



    //scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<SphereShape>(glm::vec3(0,0,0),20),nullptr,nullptr,outsideMedium));


    /*
    auto lightFunc = [](const Ray& ray){
        float a = 0.5f * (ray.dir.y + 1.0f);
        return 3.f * ((1.0f - a) * glm::vec3(1, 0.85, 0.55) + a * glm::vec3(0.45, 0.65, 1));
        };
    */

    scene->infiniteLights.push_back(std::make_shared<UniformInfiniteLight>(glm::vec3(0.45, 0.65, 1)));
    //scene->infiniteLights.push_back(std::make_shared<FunctionInfiniteLight>(lightFunc));
    //scene->infiniteLights.push_back(std::make_shared<TextureInfiniteLight>(std::make_shared<ImageTexture>("/home/markov/Documents/Coding/CPP/raytracing_in_one_weekend/kloofendal_48d_partly_cloudy_puresky.jpg",true),5));
    //scene->infiniteLights.push_back(std::make_shared<TextureInfiniteLight>(std::make_shared<FloatImageTexture>("/home/markov/Downloads/kloofendal_48d_partly_cloudy_puresky_8k.hdr"),600/255.f,[](float r){return 4 * std::sqrt(r);}));

    std::shared_ptr<LightSampler> ls = std::make_shared<PowerLightSampler>();
    scene->BuildTlas<TLAS4>();
    ls->Add(scene->GetLights());
    //ls->Add(std::make_shared<PointLight>(glm::vec3(0.3,1.5,0),glm::vec3(6)));
    ls->PreProcess(scene->BoundingBox());

    double fov = 1.7;


    glm::dvec3 lookfrom = { -1000,300,0 };
    lookfrom = { 17.3,1.2,7.2 };

    lookfrom = { 0.3,0.4,1 };//dragon

    glm::dvec3 lookat = { 0,300,0 };
    lookat = { 0,0,0 };

    //64*4

    //550
    int samples = 16;//64*16*4 -> 4 hours
    int sqrts = std::sqrt(samples);

    std::shared_ptr<Film> film = std::make_shared<Film>(glm::ivec2 { 1920,1080 }, std::make_shared<MitchellFilter>());


    auto camera = std::make_shared<Camera>(lookfrom, lookat, fov, film);
    auto sampler = std::make_shared<StratifiedSampler>(sqrts, sqrts);
    auto integrator = std::make_shared<VolPathIntegrator>(scene, camera, sampler, ls, 128);

    //camera->SetMedium(outsideMedium);
    //scene->SetMedium(outsideMedium);

    integrator->Render();
    camera->GetFilm()->WritePNG("RenderedScene");
    camera->GetFilm()->WritePPM("RenderedScene");
    camera->GetFilm()->WriteJPG("RenderedScene", 100);
    ResourceManager::get_instance().releaseTextureCache();
    ResourceManager::get_instance().releaseMeshCache();
    ResourceManager::get_instance().releaseModelCache();
}



void helmet(){
    auto scene = std::make_shared<Scene>();
    auto white = ResourceManager::get_instance().GetTexture<SolidColor>("whiteTexture", glm::vec3(.9));
    auto green = ResourceManager::get_instance().GetTexture<SolidColor>("greenTexture", glm::vec3(.2, .3, .1));
    auto light = std::make_shared<MicrofacetDiffuse>(glm::vec3(0));
    auto ch = std::make_shared<MicrofacetDiffuse>(glm::vec3 { .2,.3,.1 });
    auto glass = std::make_shared<MicrofacetDielectric>(1.5, 0.0, glm::vec3(1));
    auto checker = std::make_shared<MicrofacetDiffuse>(ResourceManager::get_instance().GetTexture<CheckerTexture>("greenTexture", white, green, glm::vec2 { 0.02 }));
    ch = std::make_shared<MicrofacetDiffuse>(std::make_shared<CheckerTexture>(white, green, glm::vec2 { 0.001,0.001 }));
    std::shared_ptr<AreaLight> area = std::make_shared<AreaLight>(std::make_shared<QuadShape>(glm::vec3(0.3, 6, 0), glm::vec3(-1, 0, 0), glm::vec3(0, 0, -1)), glm::vec3(500), false);
    auto outsideMedium = std::make_shared<HomogeneusMedium>(glm::vec3 { 0.01f, 0.9f, 0.9f }, glm::vec3 { 1.0f, 0.1f, 0.1f }, std::make_shared<HenyeyGreenstein>(0.9), 0.5f);


    //scene->Add(std::make_shared<GeometricPrimitive>(area->getShape(), light, area, nullptr));//-0.3, -1
    scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<QuadShape>(glm::vec3(-100, -0.3, -100), glm::vec3(1000, 0, 0), glm::vec3(0, 0, 1000)), ch, nullptr, nullptr));

    glm::mat4 pos = glm::mat4(1);
#define HELMET 1

#define Mig23 0
#if HELMET

    pos = glm::translate(pos, { 0,2,0 });
    pos = glm::rotate(pos, glm::radians(90.0f), glm::normalize(glm::vec3(1, 0, 0)));
    pos = glm::scale(pos, glm::vec3(2, 2, 2));
    std::shared_ptr<Primitive> transformedModel = std::make_shared<TransformedPrimitive>(ResourceManager::get_instance().CacheModel<BLAS4>("Helmet", "/home/markov/Downloads/DamagedHelmet.assbin"), pos);

#elif Mig23
    pos = glm::translate(pos, { 4,-0.3,0 });
    //pos = glm::rotate(pos,glm::radians(180.0f),glm::normalize(glm::vec3(0,1,0)));

    std::shared_ptr<Primitive> transformedModel = std::make_shared<TransformedPrimitive>(ResourceManager::get_instance().CacheModel<BLAS4>("Helmet", "/home/markov/Downloads/Models/Mig23/scene.assbin"), pos);
#else
    pos = glm::translate(pos, { 0,1 - 0.24,0 });
    pos = glm::rotate(pos, glm::radians(180.0f), glm::normalize(glm::vec3(0, 1, 0)));

    std::shared_ptr<Primitive> transformedModel = std::make_shared<TransformedPrimitive>(ResourceManager::get_instance().CacheModel<BLAS4>("Helmet", "/home/markov/Downloads/Models/Tank/scene.assbin"), pos);
#endif


    scene->Add(transformedModel);

    //scene->Add(ResourceManager::get_instance().CacheModel<BLAS4>("Medium Dragon","/home/markov/Downloads/DamagedHelmet.gltf"));//was1.2

    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(0,0.1,-1.2),0.5),std::make_shared<MicrofacetDiffuse>(glm::vec3(0.1, 0.2, 0.5)),nullptr));
    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(-1,0,-1),0.5),std::make_shared<dielectric>(1.5),nullptr));
    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(-1,0,-1),0.4),std::make_shared<dielectric>(1/1.5),nullptr));
    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(1,0,-1),0.5),std::make_shared<metal>(glm::vec3(0.8, 0.6, 0.2)),nullptr));

    //scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<SphereShape>(glm::vec3(1,0,-1),0.5),nullptr,nullptr,std::make_shared<HomogeneusMedium>(glm::vec3{0.01f, 0.9f, 0.9f},glm::vec3{1.0f, 0.1f, 0.1f},std::make_shared<HenyeyGreenstein>(0.8),25.0f,glm::vec3{1,1,1},1)));//Was 1 
        //d_list[1] = new Sphere{glm::vec3(-0.8,1,-0.5), 0.5,
        //                        new Light(glm::vec3(8, 8, 8))};
    //scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<SphereShape>(glm::vec3(0,0,0),20),nullptr,nullptr,outsideMedium));


    /*
    auto lightFunc = [](const Ray& ray){
        float a = 0.5f * (ray.dir.y + 1.0f);
        return 3.f * ((1.0f - a) * glm::vec3(1, 0.85, 0.55) + a * glm::vec3(0.45, 0.65, 1));
        };
    */

    //scene->infiniteLights.push_back(std::make_shared<UniformInfiniteLight>(glm::vec3{0,0,1}));
    //scene->infiniteLights.push_back(std::make_shared<FunctionInfiniteLight>(lightFunc));
    //scene->infiniteLights.push_back(std::make_shared<TextureInfiniteLight>(std::make_shared<ImageTexture>("/home/markov/Documents/Coding/CPP/raytracing_in_one_weekend/kloofendal_48d_partly_cloudy_puresky.jpg",true),5));
    scene->infiniteLights.push_back(std::make_shared<TextureInfiniteLight>(std::make_shared<FloatImageTexture>("/home/markov/Downloads/kloofendal_48d_partly_cloudy_puresky_8k.hdr"), 600 / 255.f, [](float r){return 40 * std::sqrt(r);}));

    std::shared_ptr<LightSampler> ls = std::make_shared<PowerLightSampler>();
    scene->BuildTlas<TLAS4>();
    ls->Add(scene->GetLights());
    //ls->Add(std::make_shared<PointLight>(glm::vec3(0.3,1.5,0),glm::vec3(6)));
    ls->PreProcess(scene->BoundingBox());

#if HELMET
    double fov = 1.5;

    glm::dvec3 lookfrom = { 4,4,7 };

    glm::dvec3 lookat = { 0,2,0 };
#elif Mig23
    double fov = 1.4;

    glm::dvec3 lookfrom = { 1.5,0.2,2 };

    glm::dvec3 lookat = { 0,0.3,0 };
#else 
    double fov = 1.4;

    glm::dvec3 lookfrom = { 6,1.2,8 };

    glm::dvec3 lookat = { 0,1.3,2 };
#endif
    //64*4
    int samples = 100;//64*16*4 -> 4 hours
    //1600 is big
    //5760
    int sqrts = std::sqrt(samples);
    //16 is 10 sec * 240 == 2400 sec which is less than our
    //
#if 0
    int frames = 300;
    for(int i = 0; i < frames;i++){
        std::shared_ptr<Film> film = std::make_shared<Film>(glm::ivec2 { 1920,1080 }, std::make_shared<MitchellFilter>());
        fov = 1.4;

        lookat = { 0,2,0.5 };
        float radius = 11;
        float angle = i / (float)frames * 2 * std::numbers::pi_v<float>;
        float x = radius * std::cos(angle);
        float z = radius * std::sin(angle);
        lookfrom = { x,3.2,z };


        auto camera = std::make_shared<Camera>(lookfrom, lookat, fov, film);
        auto sampler = std::make_shared<StratifiedSampler>(sqrts, sqrts);
        auto integrator = std::make_shared<PathIntegrator>(scene, camera, sampler, ls, 64);

        std::string index = std::to_string(i);
        std::string suffix(4, '0');
        for(int k = 4 - index.size();k < 4;k++){
            suffix[k] = index[k - 4 + index.size()];
        }
        std::cout << "Rendering: " << suffix << std::endl;
        integrator->Render();

        camera->GetFilm()->WriteJPG("Video/VideoScene" + suffix, 100);
    }
    //mkv is better?
    [[maybe_unused]] system("ffmpeg -y -framerate 30 -pattern_type glob -i 'Output/Video/VideoScene*.jpg' -c:v libx264 -pix_fmt yuv420p Output/Video/Tank.mp4");
    [[maybe_unused]] system("ffmpeg -y -framerate 30 -pattern_type glob -i 'Output/Video/VideoScene*.jpg' -c:v libx264 -crf 18 -preset veryslow -pix_fmt yuv420p Output/Video/TankHQ.mp4");
    //system("ffmpeg -y -framerate 30 -pattern_type glob -i 'Output/Video/VideoScene*.jpg' -c:v libx264 -crf 0 -preset veryslow -pix_fmt yuv420p Output/Video/TankHQ2.mp4");
    //system("ffmpeg -y -framerate 30 -pattern_type glob -i 'Output/Video/VideoScene*.jpg' -c:v libx264 -crf 0 -preset veryslow -pix_fmt yuv444p Output/Video/TankHQL.mkv");
    [[maybe_unused]] system("ffmpeg -y -framerate 60 -pattern_type glob -i 'Output/Video/VideoScene*.jpg' -c:v libx264 -pix_fmt yuv420p Output/Video/Tank60.mp4");
    [[maybe_unused]] system("ffmpeg -y -framerate 60 -pattern_type glob -i 'Output/Video/VideoScene*.jpg' -c:v libx264 -crf 18 -preset veryslow -pix_fmt yuv420p Output/Video/TankHQ60.mp4");
    //system("ffmpeg -y -framerate 60 -pattern_type glob -i 'Output/Video/VideoScene*.jpg' -c:v libx264 -crf 0 -preset veryslow -pix_fmt yuv420p Output/Video/TankHQ260.mp4");
    //system("ffmpeg -y -framerate 60 -pattern_type glob -i 'Output/Video/VideoScene*.jpg' -c:v libx264 -crf 0 -preset veryslow -pix_fmt yuv444p Output/Video/TankHQL60.mkv");
    //system("rm Output/Video/VideoScene*.jpg");
    //camera->GetFilm()->WritePNG("RenderedScene");
    //camera->GetFilm()->WritePPM("RenderedScene");
#else
    std::shared_ptr<Film> film = std::make_shared<Film>(glm::ivec2 { 1920,1080 }, std::make_shared<MitchellFilter>());




    auto camera = std::make_shared<Camera>(lookfrom, lookat, fov, film);
    auto sampler = std::make_shared<StratifiedSampler>(sqrts, sqrts);
    auto integrator = std::make_shared<PathIntegrator>(scene, camera, sampler, ls, 64);

    integrator->Render();

    camera->GetFilm()->WriteJPG("Helmet", 100);
#endif
    ResourceManager::get_instance().releaseTextureCache();
    ResourceManager::get_instance().releaseMeshCache();
    ResourceManager::get_instance().releaseModelCache();
}



void knight(){
    auto scene = std::make_shared<Scene>();
    auto white = ResourceManager::get_instance().GetTexture<SolidColor>("whiteTexture", glm::vec3(.9));
    auto green = ResourceManager::get_instance().GetTexture<SolidColor>("greenTexture", glm::vec3(.2, .3, .1));
    auto light = std::make_shared<MicrofacetDiffuse>(glm::vec3(0));
    auto ch = std::make_shared<MicrofacetDiffuse>(glm::vec3 { .2,.3,.1 });
    auto glass = std::make_shared<MicrofacetDielectric>(1.5, 0.0, glm::vec3(1));
    auto checker = std::make_shared<MicrofacetDiffuse>(ResourceManager::get_instance().GetTexture<CheckerTexture>("greenTexture", white, green, glm::vec2 { 0.02 }));
    ch = std::make_shared<MicrofacetDiffuse>(std::make_shared<CheckerTexture>(white, green, glm::vec2 { 0.001,0.001 }));
    std::shared_ptr<AreaLight> area = std::make_shared<AreaLight>(std::make_shared<QuadShape>(glm::vec3(0.3, 6, 0), glm::vec3(-1, 0, 0), glm::vec3(0, 0, -1)), glm::vec3(500), false);
    auto outsideMedium = std::make_shared<HomogeneusMedium>(glm::vec3 { 0.01f, 0.9f, 0.9f }, glm::vec3 { 1.0f, 0.1f, 0.1f }, std::make_shared<HenyeyGreenstein>(0.9), 0.5f);


    //scene->Add(std::make_shared<GeometricPrimitive>(area->getShape(), light, area, nullptr));//-0.3, -1
    //scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<QuadShape>(glm::vec3(-100,-0.3,-100), glm::vec3(1000,0,0), glm::vec3(0,0,1000)), ch, nullptr, nullptr));

    glm::mat4 pos = glm::mat4(1);

    pos = glm::translate(pos, { 0,0,0 });
    pos = glm::rotate(pos, glm::radians(180.0f), glm::normalize(glm::vec3(0, 1, 0)));

    std::shared_ptr<Primitive> transformedModel = std::make_shared<TransformedPrimitive>(ResourceManager::get_instance().CacheModel<BLAS4>("Helmet", "/home/markov/Downloads/AttenuationTest.gltf"), pos);
    //DragonAttenuation



    scene->Add(transformedModel);

    //scene->Add(ResourceManager::get_instance().CacheModel<BLAS4>("Medium Dragon","/home/markov/Downloads/DamagedHelmet.gltf"));//was1.2

    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(0,0.1,-1.2),0.5),std::make_shared<MicrofacetDiffuse>(glm::vec3(0.1, 0.2, 0.5)),nullptr));
    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(-1,0,-1),0.5),std::make_shared<dielectric>(1.5),nullptr));
    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(-1,0,-1),0.4),std::make_shared<dielectric>(1/1.5),nullptr));
    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(1,0,-1),0.5),std::make_shared<metal>(glm::vec3(0.8, 0.6, 0.2)),nullptr));

    //scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<SphereShape>(glm::vec3(1,0,-1),0.5),nullptr,nullptr,std::make_shared<HomogeneusMedium>(glm::vec3{0.01f, 0.9f, 0.9f},glm::vec3{1.0f, 0.1f, 0.1f},std::make_shared<HenyeyGreenstein>(0.8),25.0f,glm::vec3{1,1,1},1)));//Was 1 
        //d_list[1] = new Sphere{glm::vec3(-0.8,1,-0.5), 0.5,
        //                        new Light(glm::vec3(8, 8, 8))};
    //scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<SphereShape>(glm::vec3(0,0,0),20),nullptr,nullptr,outsideMedium));

    /*
    auto lightFunc = [](const Ray& ray){
        float a = 0.5f * (ray.dir.y + 1.0f);
        return 3.f * ((1.0f - a) * glm::vec3(1, 0.85, 0.55) + a * glm::vec3(0.45, 0.65, 1));
        };
    */

    //scene->infiniteLights.push_back(std::make_shared<UniformInfiniteLight>(glm::vec3{0,0,1}));
    //scene->infiniteLights.push_back(std::make_shared<FunctionInfiniteLight>(lightFunc));
    //scene->infiniteLights.push_back(std::make_shared<TextureInfiniteLight>(std::make_shared<ImageTexture>("/home/markov/Documents/Coding/CPP/raytracing_in_one_weekend/kloofendal_48d_partly_cloudy_puresky.jpg",true),5));
    //scene->infiniteLights.push_back(std::make_shared<TextureInfiniteLight>(std::make_shared<FloatImageTexture>("/home/markov/Downloads/kloofendal_48d_partly_cloudy_puresky_8k.hdr"),600/255.f,[](float r){return 40 * std::sqrt(r);}));
    scene->infiniteLights.push_back(std::make_shared<TextureInfiniteLight>(std::make_shared<FloatImageTexture>("/home/markov/Downloads/artist_workshop_2k.hdr"), 1, [](float r){return 40 * std::sqrt(r);}));


    std::shared_ptr<LightSampler> ls = std::make_shared<PowerLightSampler>();
    scene->BuildTlas<TLAS4>();
    ls->Add(scene->GetLights());
    //ls->Add(std::make_shared<PointLight>(glm::vec3(0.3,1.5,0),glm::vec3(6)));
    ls->PreProcess(scene->BoundingBox());


    double fov = 0.45;
    fov = 1.1;
    glm::dvec3 lookfrom = { 0,1,10 };
    lookfrom = { 0,0,-15 };//lookfrom = 0,1,0 NaN ? 
    glm::dvec3 lookat = { 0,0.5,0 };
    lookat = { 0,0,0 };

    int samples = 100 * 4;
    int sqrts = std::sqrt(samples);

    //16 is 10 sec * 240 == 2400 sec which is less than our
    //

    std::shared_ptr<Film> film = std::make_shared<Film>(glm::ivec2 { 1000,1000 }, std::make_shared<MitchellFilter>());


    auto camera = std::make_shared<Camera>(lookfrom, lookat, fov, film);
    auto sampler = std::make_shared<StratifiedSampler>(sqrts, sqrts);
    auto integrator = std::make_shared<VolPathIntegrator>(scene, camera, sampler, ls, 64);

    integrator->Render();

    camera->GetFilm()->WriteJPG("Knight", 100);


    //camera->GetFilm()->WritePNG("RenderedScene");
    //camera->GetFilm()->WritePPM("RenderedScene");
    ResourceManager::get_instance().releaseTextureCache();
    ResourceManager::get_instance().releaseMeshCache();
    ResourceManager::get_instance().releaseModelCache();
}

void opacity(){
    auto scene = std::make_shared<Scene>();
    auto white = ResourceManager::get_instance().GetTexture<SolidColor>("whiteTexture", glm::vec3(.9));
    auto green = ResourceManager::get_instance().GetTexture<SolidColor>("greenTexture", glm::vec3(.2, .3, .1));
    auto light = std::make_shared<MicrofacetDiffuse>(glm::vec3(0));
    auto ch = std::make_shared<MicrofacetDiffuse>(glm::vec3 { .2,.3,.1 });
    auto glass = std::make_shared<MicrofacetDielectric>(1.5, 0.0, glm::vec3(1));
    auto checker = std::make_shared<MicrofacetDiffuse>(ResourceManager::get_instance().GetTexture<CheckerTexture>("greenTexture", white, green, glm::vec2 { 0.02 }));
    ch = std::make_shared<MicrofacetDiffuse>(std::make_shared<CheckerTexture>(white, green, glm::vec2 { 0.001,0.001 }));
    std::shared_ptr<AreaLight> area = std::make_shared<AreaLight>(std::make_shared<QuadShape>(glm::vec3(0.3, 6, 0), glm::vec3(-1, 0, 0), glm::vec3(0, 0, -1)), glm::vec3(500), false);
    auto outsideMedium = std::make_shared<HomogeneusMedium>(glm::vec3 { 0.01f, 0.9f, 0.9f }, glm::vec3 { 1.0f, 0.1f, 0.1f }, std::make_shared<HenyeyGreenstein>(0.9), 0.5f);


    //scene->Add(std::make_shared<GeometricPrimitive>(area->getShape(), light, area, nullptr));//-0.3, -1
    //scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<QuadShape>(glm::vec3(-100,-0.3,-100), glm::vec3(1000,0,0), glm::vec3(0,0,1000)), ch, nullptr, nullptr));

    glm::mat4 pos = glm::mat4(1);

    pos = glm::translate(pos, { 0,0,0 });
    pos = glm::rotate(pos, glm::radians(180.0f), glm::normalize(glm::vec3(0, 1, 0)));

    std::shared_ptr<Primitive> transformedModel = std::make_shared<TransformedPrimitive>(ResourceManager::get_instance().CacheModel<BLAS4>("Helmet", "/home/markov/Downloads/TransmissionRoughnessTest.gltf"), pos);
    //DragonAttenuation



    scene->Add(transformedModel);

    //scene->Add(ResourceManager::get_instance().CacheModel<BLAS4>("Medium Dragon","/home/markov/Downloads/DamagedHelmet.gltf"));//was1.2

    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(0,0.1,-1.2),0.5),std::make_shared<MicrofacetDiffuse>(glm::vec3(0.1, 0.2, 0.5)),nullptr));
    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(-1,0,-1),0.5),std::make_shared<dielectric>(1.5),nullptr));
    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(-1,0,-1),0.4),std::make_shared<dielectric>(1/1.5),nullptr));
    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(1,0,-1),0.5),std::make_shared<metal>(glm::vec3(0.8, 0.6, 0.2)),nullptr));

    //scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<SphereShape>(glm::vec3(1,0,-1),0.5),nullptr,nullptr,std::make_shared<HomogeneusMedium>(glm::vec3{0.01f, 0.9f, 0.9f},glm::vec3{1.0f, 0.1f, 0.1f},std::make_shared<HenyeyGreenstein>(0.8),25.0f,glm::vec3{1,1,1},1)));//Was 1 
        //d_list[1] = new Sphere{glm::vec3(-0.8,1,-0.5), 0.5,
        //                        new Light(glm::vec3(8, 8, 8))};
    //scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<SphereShape>(glm::vec3(0,0,0),20),nullptr,nullptr,outsideMedium));

    auto lightFunc = [](const Ray& ray){
        float a = 0.5f * (ray.dir.y + 1.0f);
        return 3.f * ((1.0f - a) * glm::vec3(1, 0.85, 0.55) + a * glm::vec3(0.45, 0.65, 1));
        };

    //scene->infiniteLights.push_back(std::make_shared<UniformInfiniteLight>(glm::vec3{0,0,1}));
    scene->infiniteLights.push_back(std::make_shared<FunctionInfiniteLight>(lightFunc));
    //scene->infiniteLights.push_back(std::make_shared<TextureInfiniteLight>(std::make_shared<ImageTexture>("/home/markov/Documents/Coding/CPP/raytracing_in_one_weekend/kloofendal_48d_partly_cloudy_puresky.jpg",true),5));
    //scene->infiniteLights.push_back(std::make_shared<TextureInfiniteLight>(std::make_shared<FloatImageTexture>("/home/markov/Downloads/kloofendal_48d_partly_cloudy_puresky_8k.hdr"),600/255.f,[](float r){return 40 * std::sqrt(r);}));
    //scene->infiniteLights.push_back(std::make_shared<TextureInfiniteLight>(std::make_shared<FloatImageTexture>("/home/markov/Downloads/artist_workshop_2k.hdr"),1,[](float r){return 40 * std::sqrt(r);}));


    std::shared_ptr<LightSampler> ls = std::make_shared<PowerLightSampler>();
    scene->BuildTlas<TLAS4>();
    ls->Add(scene->GetLights());
    //ls->Add(std::make_shared<PointLight>(glm::vec3(0.3,1.5,0),glm::vec3(6)));
    ls->PreProcess(scene->BoundingBox());


    double fov = 0.45;
    fov = 1.3;
    glm::dvec3 lookfrom = { 0,1,10 };
    lookfrom = { 0.05,-0.05,-1 };//lookfrom = 0,1,0 NaN ? 
    glm::dvec3 lookat = { 0,0.5,0 };
    lookat = { 0.05,-0.05,0 };

    int samples = 100;
    int sqrts = std::sqrt(samples);

    //16 is 10 sec * 240 == 2400 sec which is less than our
    //

    std::shared_ptr<Film> film = std::make_shared<Film>(glm::ivec2 { 1920,1080 }, std::make_shared<MitchellFilter>());


    auto camera = std::make_shared<Camera>(lookfrom, lookat, fov, film);
    auto sampler = std::make_shared<StratifiedSampler>(sqrts, sqrts);
    auto integrator = std::make_shared<VolPathIntegrator>(scene, camera, sampler, ls, 64);

    integrator->Render();

    camera->GetFilm()->WriteJPG("Knight", 100);


    //camera->GetFilm()->WritePNG("RenderedScene");
    //camera->GetFilm()->WritePPM("RenderedScene");
    ResourceManager::get_instance().releaseTextureCache();
    ResourceManager::get_instance().releaseMeshCache();
    ResourceManager::get_instance().releaseModelCache();
}

void transmission(){
    auto scene = std::make_shared<Scene>();
    auto white = ResourceManager::get_instance().GetTexture<SolidColor>("whiteTexture", glm::vec3(.9));
    auto green = ResourceManager::get_instance().GetTexture<SolidColor>("greenTexture", glm::vec3(.2, .3, .1));
    auto light = std::make_shared<MicrofacetDiffuse>(glm::vec3(0));
    auto ch = std::make_shared<MicrofacetDiffuse>(glm::vec3 { .2,.3,.1 });
    auto glass = std::make_shared<MicrofacetDielectric>(1.5, 0.0, glm::vec3(1));
    auto checker = std::make_shared<MicrofacetDiffuse>(ResourceManager::get_instance().GetTexture<CheckerTexture>("greenTexture", white, green, glm::vec2 { 0.02 }));
    ch = std::make_shared<MicrofacetDiffuse>(std::make_shared<CheckerTexture>(white, green, glm::vec2 { 0.001,0.001 }));
    std::shared_ptr<AreaLight> area = std::make_shared<AreaLight>(std::make_shared<QuadShape>(glm::vec3(0.3, 6, 0), glm::vec3(-1, 0, 0), glm::vec3(0, 0, -1)), glm::vec3(500), false);
    auto outsideMedium = std::make_shared<HomogeneusMedium>(glm::vec3 { 0.01f, 0.9f, 0.9f }, glm::vec3 { 1.0f, 0.1f, 0.1f }, std::make_shared<HenyeyGreenstein>(0.9), 0.5f);


    //scene->Add(std::make_shared<GeometricPrimitive>(area->getShape(), light, area, nullptr));//-0.3, -1
    //scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<QuadShape>(glm::vec3(-100,-0.3,-100), glm::vec3(1000,0,0), glm::vec3(0,0,1000)), ch, nullptr, nullptr));

    glm::mat4 pos = glm::mat4(1);

    pos = glm::translate(pos, { 0,0,0 });
    pos = glm::rotate(pos, glm::radians(180.0f), glm::normalize(glm::vec3(0, 1, 0)));

    std::shared_ptr<Primitive> transformedModel = std::make_shared<TransformedPrimitive>(ResourceManager::get_instance().CacheModel<BLAS4>("Helmet", "/home/markov/Downloads/AlphaBlendModeTest.gltf"), pos);
    //DragonAttenuation



    scene->Add(transformedModel);

    //scene->Add(ResourceManager::get_instance().CacheModel<BLAS4>("Medium Dragon","/home/markov/Downloads/DamagedHelmet.gltf"));//was1.2

    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(0,0.1,-1.2),0.5),std::make_shared<MicrofacetDiffuse>(glm::vec3(0.1, 0.2, 0.5)),nullptr));
    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(-1,0,-1),0.5),std::make_shared<dielectric>(1.5),nullptr));
    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(-1,0,-1),0.4),std::make_shared<dielectric>(1/1.5),nullptr));
    //scene->Add(new GeometricPrimitive(new SphereShape(glm::vec3(1,0,-1),0.5),std::make_shared<metal>(glm::vec3(0.8, 0.6, 0.2)),nullptr));

    //scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<SphereShape>(glm::vec3(1,0,-1),0.5),nullptr,nullptr,std::make_shared<HomogeneusMedium>(glm::vec3{0.01f, 0.9f, 0.9f},glm::vec3{1.0f, 0.1f, 0.1f},std::make_shared<HenyeyGreenstein>(0.8),25.0f,glm::vec3{1,1,1},1)));//Was 1 
        //d_list[1] = new Sphere{glm::vec3(-0.8,1,-0.5), 0.5,
        //                        new Light(glm::vec3(8, 8, 8))};
    //scene->Add(std::make_shared<GeometricPrimitive>(std::make_shared<SphereShape>(glm::vec3(0,0,0),20),nullptr,nullptr,outsideMedium));

    auto lightFunc = [](const Ray& ray){
        float a = 0.5f * (ray.dir.y + 1.0f);
        return 3.f * ((1.0f - a) * glm::vec3(1, 0.85, 0.55) + a * glm::vec3(0.45, 0.65, 1));
        };

    //scene->infiniteLights.push_back(std::make_shared<UniformInfiniteLight>(glm::vec3{0,0,1}));
    scene->infiniteLights.push_back(std::make_shared<FunctionInfiniteLight>(lightFunc));
    //scene->infiniteLights.push_back(std::make_shared<TextureInfiniteLight>(std::make_shared<ImageTexture>("/home/markov/Documents/Coding/CPP/raytracing_in_one_weekend/kloofendal_48d_partly_cloudy_puresky.jpg",true),5));
    //scene->infiniteLights.push_back(std::make_shared<TextureInfiniteLight>(std::make_shared<FloatImageTexture>("/home/markov/Downloads/kloofendal_48d_partly_cloudy_puresky_8k.hdr"),600/255.f,[](float r){return 40 * std::sqrt(r);}));
    //scene->infiniteLights.push_back(std::make_shared<TextureInfiniteLight>(std::make_shared<FloatImageTexture>("/home/markov/Downloads/artist_workshop_2k.hdr"),1,[](float r){return 40 * std::sqrt(r);}));


    std::shared_ptr<LightSampler> ls = std::make_shared<PowerLightSampler>();
    scene->BuildTlas<TLAS4>();
    ls->Add(scene->GetLights());
    //ls->Add(std::make_shared<PointLight>(glm::vec3(0.3,1.5,0),glm::vec3(6)));
    ls->PreProcess(scene->BoundingBox());


    double fov = 0.45;
    fov = 1.3;
    glm::dvec3 lookfrom = { 0,1,10 };
    lookfrom = { 0.05,-0.05,-1 };//lookfrom = 0,1,0 NaN ? 
    lookfrom = { 0,0,-7 };
    glm::dvec3 lookat = { 0,0.5,0 };
    lookat = { 0.05,-0.05,0 };

    int samples = 36;
    int sqrts = std::sqrt(samples);

    //16 is 10 sec * 240 == 2400 sec which is less than our
    //

    std::shared_ptr<Film> film = std::make_shared<Film>(glm::ivec2 { 1920,1080 }, std::make_shared<MitchellFilter>());


    auto camera = std::make_shared<Camera>(lookfrom, lookat, fov, film);
    auto sampler = std::make_shared<StratifiedSampler>(sqrts, sqrts);
    auto integrator = std::make_shared<VolPathIntegrator>(scene, camera, sampler, ls, 64);

    integrator->Render();

    camera->GetFilm()->WriteJPG("Knight", 100);


    //camera->GetFilm()->WritePNG("RenderedScene");
    //camera->GetFilm()->WritePPM("RenderedScene");
    ResourceManager::get_instance().releaseTextureCache();
    ResourceManager::get_instance().releaseMeshCache();
    ResourceManager::get_instance().releaseModelCache();
}
//how to deal when we spawn ray inside objects with medium?
//interaction spawnRay() -> this spawns ray but gives us 
// ### />   <- .
//getMedium would be ###, but we would set medium to be null (coming from outside)
// </###   <- .
//getMedium would be nullptr, but we would set medium to be ### (coming from inside)

//examples folder -> give links to models

int main(){
    //stbi_set_flip_vertically_on_load(true);
    //"/home/markov/Documents/Coding/CPP/testing/stanford/common-3d-test-models-master/data/lucy.obj"
    switch(1){
    case 0:
        temp();
        break;
    case 1:
        Miguel();
        break;
    case 2:
        NoModel();
        break;
    case 3:
        MatTest();
        break;
    case 4:
        helmet();
        break;
    case 5:
        knight();
        break;
    case 6:
        opacity();
        break;
    case 7:
        transmission();
        break;
    }
    ResourceManager::get_instance().releaseTextureCache();
    ResourceManager::get_instance().releaseMeshCache();
    ResourceManager::get_instance().releaseModelCache();
}
