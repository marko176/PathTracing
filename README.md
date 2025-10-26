A physically based offline CPU pathtracer written in C++20

 - [Screenshots](#screenshots)
 - [Requirements](#requirements)
 - [Building](#build)
 - [Features](#features)
 - [Dependencies](#dependencies)
   
## Screenshots

San Miguel - Rendered with 1024 SPP, 128 bounces
![san_miguel](assets/SanMiguel1024.jpg)

## REQUIREMENTS

- **C++20** compiler
- **Cmake 3.22+** for building

## Build

1. Cloning the project with submodules:
```bash
git clone --recurse-submodules https://github.com/marko176/PathTracing.git
```
2. Building the project:
```bash
cd PathTracing
mkdir build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build ./build --config Release -j N   (set N to thread count)
```
3. Running the pathtracer do: `./Release/PathTracer`

Output images will be stored in: `Output/`

## Features:

- **Camera** with field of view and thin-lens approximation

- **BVH2** using SAH-based binned building and simd during traversal, split into:
  - TLAS (Top-Level Acceleration Structure)
  - BLAS (Bottom-Level Acceleration Structure)

- **Model Loading** with [Assimp](https://github.com/assimp/assimp)
  - Supports obj and gltf formats
  - Supports model serialization for faster reloads in assbin format

- **Physically Based Materials**
  - Microfacet Dielectric
  - Microfacet Diffuse
  - Specular Dielectric
  - Specular Conductor

- **Volumetric Rendering**
  - Supports **participating media** (homogeneous)  
  - Handles volumetric scattering and absorption within the path tracer  
  - Works with MIS and NEE for efficient light transport through media

- **Sampling**
  - Cosine-weighted importance sampling for Lambertian BRDF
  - [VNDF](https://jcgt.org/published/0007/04/01/) (visible normal distribution function) importance sampling for GGX specular BRDF
  - **Multiple Importance Sampling (MIS)** with Next Event Estimation (NEE) weighted with the power heuristic [Veach (1997)](https://graphics.stanford.edu/papers/veach_thesis/)
  - **Russian Roulette** path termination

- **Texture Mapping**
  - Albedo, Normal, Roughness, Metallic, Alpha, Emissive
  - support for HDR textures

## Dependencies
  - [Assimp](https://github.com/assimp/assimp)
  - [stb_image](https://github.com/nothings/stb)
  - [glm](https://github.com/g-truc/glm)
  - [pcg](https://github.com/imneme/pcg-cpp)



