A simple C++ pathtracer

**REQUIREMENTS**:

- **C++20** compiler
- **Cmake 3.22+** for building


Cloning the project:
```bash
git clone --recurse-submodules https://github.com/marko176/PathTracing.git
```
Building the project:
```bash
mkdir build
cmake -S . -B build
cmake --build ./build -j N   (set N to thread count)
```
To run the pathtracer do: `./PathTracer.exe`
