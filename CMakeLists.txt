cmake_minimum_required(VERSION 3.22)
project(PathTracer)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")



set(CMAKE_CXX_COMPILER "/usr/bin/x86_64-linux-gnu-g++-14")
set(CMAKE_CXX_STANDARD 23) 
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O3 -g -march=native")
set(INTERPROCEDURAL_OPTIMIZATION TRUE)
set(BUILD_SHARED_LIBS OFF)
add_subdirectory(assimp)

add_executable(PathTracer 
            stb_image.cpp 
            main_BVH.cpp 
            Mesh.cpp 
            Model.cpp 
            Texture.cpp 
            ResourceManager.cpp 
            Primitive.cpp 
            Shape.cpp
            Light.cpp
            LightSampler.cpp)

target_link_libraries(PathTracer PRIVATE assimp)
target_include_directories(PathTracer PRIVATE ${CMAKE_SOURCE_DIR})
