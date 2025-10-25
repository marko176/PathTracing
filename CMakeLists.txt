cmake_minimum_required(VERSION 3.22)
project(PathTracer)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/$<CONFIG>")


set(BUILD_SHARED_LIBS OFF)
set(ASSIMP_BUILD_TESTS OFF)

add_subdirectory(assimp)

#set(CMAKE_CXX_COMPILER "/usr/bin/x86_64-linux-gnu-g++-14")
set(CMAKE_CXX_STANDARD 20) 
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()



set(INTERPROCEDURAL_OPTIMIZATION TRUE)

add_executable(PathTracer 
            stb_image.cpp 
            stb_image_write.cpp
            main.cpp  
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


target_include_directories(PathTracer PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(PathTracer PRIVATE assimp)

if (MSVC) 
    #add_compile_options(/FS /O2)
    target_compile_options(PathTracer PRIVATE
        $<$<CONFIG:Debug>:/FS /Zi /Od /W3 /MDd>
        $<$<CONFIG:RelWithDebInfo>:/FS /Zi /O2 /W4 /MD>
        $<$<CONFIG:Release>:/FS /O2 /Ox /W4 /MD>
    )
else()
    #add_compile_options(-O3 -march=native -g)
    target_compile_options(PathTracer PRIVATE
        # Debug: no optimization, include debug info, be very noisy with warnings
        $<$<CONFIG:Debug>:-O0;-g;-march=native;-Wall;-Wshadow>
        # RelWithDebInfo: useful balance (optimized + debug info)
        $<$<CONFIG:RelWithDebInfo>:-O2;-g;-march=native;-Wall>
        # Release: max optimization, strip debug, enable LTO
        $<$<CONFIG:Release>:-O3;-g;-march=native>
    )
endif()
