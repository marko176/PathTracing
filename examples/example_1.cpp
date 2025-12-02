#include <optional>
#include <memory>
#include "../Model.hpp"
#include "../Material.hpp"
#include "../Primitive.hpp"
#include "../Light.hpp"
#include "../LightSampler.hpp"
#include "../Sampler.hpp"
#include "../Filter.hpp"
#include "../Film.hpp"
#include "../Camera.hpp"
#include "../Scene.hpp"
#include "../ResourceManager.hpp"
#include "../Integrators.hpp"


int main(){
    // 1) Scene
    // The Scene object holds all the objects (primitives) in the scene we will render
    auto scene = std::make_shared<Scene>();


    // 2) Textures
    // Create/fetch textures from the ResourceManager. The names are used as keys.
    std::shared_ptr<Texture> white = ResourceManager::get_instance().GetTexture<SolidColor>("whiteTexture", glm::vec3(.9));
    std::shared_ptr<Texture> green = ResourceManager::get_instance().GetTexture<SolidColor>("greenTexture", glm::vec3(.2, .3, .1));


    // 3) Materials
    // Materials are made from textures
    std::shared_ptr<Material> sheckerMaterial = std::make_shared<MicrofacetDiffuse>(std::make_shared<CheckerTexture>(white, green, glm::vec2 { 0.001,0.001 }));
    std::shared_ptr<Material> sphereMaterial = std::make_shared<MicrofacetDiffuse>(green);


    // 4) Shapes
    // Create a large quad for the floor and a sphere 
    std::shared_ptr<Shape> floor = std::make_shared<QuadShape>(glm::vec3(-100, -0.3, -100), glm::vec3(1000, 0, 0), glm::vec3(0, 0, 1000));
    std::shared_ptr<Shape> sphere = std::make_shared<SphereShape>(glm::vec3(0, 0.1, -1.2), 0.5);
    std::shared_ptr<Shape> mediumSphere = std::make_shared<SphereShape>(glm::vec3(1, 0, -1), 0.5);
    std::shared_ptr<Shape> lightShape = std::make_shared<QuadShape>(glm::vec3(-1, -0.28, -1), glm::vec3(0.2, 0, -0.2), glm::vec3(0, 0.2, 0));

    glm::vec3 lightColor = glm::vec3(1, 0, 0);
    float lightIntenisty = 600;
    bool oneSided = false;
    std::shared_ptr<AreaLight> areaLight = std::make_shared<AreaLight>(lightShape, lightColor * lightIntenisty, oneSided);


    // 5) Primitives
    // Primitves wrap shapes and give them additional info (materials)
    std::shared_ptr<GeometricPrimitive> floorPrimitive = std::make_shared<GeometricPrimitive>(floor, sheckerMaterial, nullptr, nullptr);
    std::shared_ptr<GeometricPrimitive> spherePrimitive = std::make_shared<GeometricPrimitive>(sphere, sphereMaterial, nullptr);
    std::shared_ptr<GeometricPrimitive> lightPrimitive = std::make_shared<GeometricPrimitive>(areaLight->getShape(), std::make_shared<MicrofacetDiffuse>(glm::vec3(0)), areaLight, nullptr);


    // Medium creation
    std::shared_ptr<PhaseFunction> phase = std::make_shared<HenyeyGreenstein>(0.8);
    float density = 5.0f;
    std::shared_ptr<Medium> homogeneusMedium = std::make_shared<HomogeneusMedium>(glm::vec3 { 0.01f, 0.9f, 0.9f }, glm::vec3 { 1.0f, 0.1f, 0.1f }, phase, density);
    std::shared_ptr<GeometricPrimitive> mediumPrimitive = std::make_shared<GeometricPrimitive>(mediumSphere, nullptr, nullptr, homogeneusMedium);


    // 6) Add primitives to the sceen
    scene->Add(floorPrimitive);
    scene->Add(spherePrimitive);
    scene->Add(lightPrimitive);
    scene->Add(mediumPrimitive);


    // 7) Create background color
    scene->infiniteLights.push_back(std::make_shared<UniformInfiniteLight>(glm::vec3(0.45, 0.65, 1)));


    // 8) Preprocess the scene
    // BuildTlas has to take some TLAS eg TLAS,TLAS4,TLAS8
    scene->BuildTlas<TLAS4>();


    // 9) Create your Film (image buffer) + filter
    auto film = std::make_shared<Film>(glm::ivec2 { 1920,1080 }, std::make_shared<MitchellFilter>());


    // 10) Camera setup
    float fov = 1.7;

    glm::dvec3 lookfrom = { 0.3,0.4,1 };

    glm::dvec3 lookat = { 0,0,0 };

    auto camera = std::make_shared<Camera>(lookfrom, lookat, fov, film);


    // 11) Create Sampler and Integrator
    int samples = 100;
    auto sampler = std::make_shared<UniformSampler>(samples);
    auto simpleIntegrator = std::make_shared<SimplePathIntegrator>(scene, camera, sampler, 128);


    // 12) Create light sampler for more advanced integrators
    auto lightSampler = std::make_shared<UniformLightSampler>();
    lightSampler->Add(scene->GetLights());
    lightSampler->PreProcess(scene->BoundingBox());
    auto pathIntegrator = std::make_shared<PathIntegrator>(scene, camera, sampler, lightSampler, 128);
    auto volumetricPathIntegrator = std::make_shared<VolPathIntegrator>(scene, camera, sampler, lightSampler, 128);


    // 13) Render the image
    simpleIntegrator->Render();
    camera->GetFilm()->WriteJPG("example_1_simple", 100);
    film->Clear();
    pathIntegrator->Render();
    camera->GetFilm()->WriteJPG("example_1_path", 100);
    film->Clear();
    volumetricPathIntegrator->Render();
    camera->GetFilm()->WriteJPG("example_1_volPath", 100);
    film->Clear();


    // Release cached resources (ResourceManager destructor also does this)
    ResourceManager::get_instance().releaseTextureCache();
    ResourceManager::get_instance().releaseMeshCache();
    ResourceManager::get_instance().releaseModelCache();
}