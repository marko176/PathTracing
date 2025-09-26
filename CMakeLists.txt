cmake_minimum_required(VERSION 3.22)
project(PathTracer)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/$<CONFIG>")


set(BUILD_SHARED_LIBS OFF)
set(ASSIMP_BUILD_TESTS OFF)

add_subdirectory(assimp)

#set(CMAKE_CXX_COMPILER "/usr/bin/x86_64-linux-gnu-g++-14")
set(CMAKE_CXX_STANDARD 23) 
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")

if (MSVC) 
    add_compile_options(/FS /O2)
else()
    add_compile_options(-O3 -march=native -g)
endif()

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
            Integrators.cpp
            PhaseFunction.cpp)

target_include_directories(PathTracer.exe PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(PathTracer.exe PRIVATE assimp)
