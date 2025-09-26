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
cd PathTracing
mkdir build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build ./build --config Release -j N   (set N to thread count)
```
To run the pathtracer do: `./Release/PathTracer`

Output images will be stored in: `Output/`