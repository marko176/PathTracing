cmake_minimum_required(VERSION 3.22)
project(PathTracer)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")


set(BUILD_SHARED_LIBS OFF)
set(ASSIMP_BUILD_TESTS OFF)

add_subdirectory(assimp)

set(CMAKE_CXX_COMPILER "/usr/bin/x86_64-linux-gnu-g++-14")
set(CMAKE_CXX_STANDARD 23) 
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O3 -g -march=native ")

set(INTERPROCEDURAL_OPTIMIZATION TRUE)

add_executable(PathTracer.exe 
            stb_image.cpp 
            main.cpp 
            Mesh.cpp 
            Model.cpp 
            Texture.cpp 
            ResourceManager.cpp 
            Primitive.cpp 
            Shape.cpp
            Light.cpp
            LightSampler.cpp
            Scene.cpp
            Integrators.cpp)

target_include_directories(PathTracer.exe PRIVATE ${CMAKE_SOURCE_DIR} "/usr/local/include/OpenImageDenoise")
target_link_libraries(PathTracer.exe PRIVATE assimp "/usr/local/lib/libOpenImageDenoise.so")
