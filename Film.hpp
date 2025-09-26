#pragma once
#include "Filter.hpp"
#include <atomic>
#include <fstream>
#include <numbers>
inline double luminance(const glm::dvec3& v){
    return dot(v, glm::dvec3(0.2126f, 0.7152f, 0.0722f));
}

inline glm::dvec3 change_luminance(const glm::dvec3& c_in, double l_out){
    double l_in = luminance(c_in);
    return c_in * (l_out / l_in);
}

inline glm::dvec3 reinhard_extended_luminance(const glm::dvec3& v, double max_white_l){
    double l_old = luminance(v);
    double numerator = l_old * (1.0 + (l_old / (max_white_l * max_white_l)));
    double l_new = numerator / (1.0 + l_old);
    return change_luminance(v, l_new);
}

inline glm::dvec3 reinhard_jodie(const glm::dvec3& v){
    double l = luminance(v);
    glm::dvec3 tv = v / (1.0 + v);
    return glm::mix(v / (1.0 + l), tv, tv);
}

inline glm::dvec3 ACESFilm(glm::dvec3 color) {
    const double A = 2.51f;
    const double B = 0.03f;
    const double C = 2.43f;
    const double D = 0.59f;
    const double E = 0.14f;
    return glm::clamp((color * (A * color + B)) / (color * (C * color + D) + E), 0.0, 1.0);
}

class Film;

struct FilmTilePixel {
    glm::dvec3 RGB = {0,0,0};
    double weight = 0;
};

class FilmTile {
public: 

    FilmTile(const Bounds2i& pixelBounds,const Bounds2i& tileBounds, const std::shared_ptr<Filter>& filter) : pixelBounds(pixelBounds), tileBounds(tileBounds) , filter(filter), radius(glm::ceil(filter->Radius() - glm::vec2{0.5,0.5})), oneOverFilterIntegral(1.0 / filter->Integral()){
        std::size_t tileWidth = tileBounds.max.x - tileBounds.min.x;
        std::size_t tileHeight = tileBounds.max.y - tileBounds.min.y;
        pixels = std::vector<FilmTilePixel>(tileWidth * tileHeight);
    }

    void Add(const glm::dvec2& p, const glm::dvec3& RGB) {

        glm::dvec2 pixelSample = glm::fract(p);
        glm::ivec2 pixel = glm::floor(p);

        std::size_t tileWidth = tileBounds.max.x - tileBounds.min.x;
        for(int y = -radius.y;y<=radius.y;y++){
            for(int x = -radius.x;x<=radius.x;x++){
                glm::dvec2 sample_pos = glm::dvec2{x,y} + glm::dvec2{0.5,0.5} - pixelSample;
                double weight = filter->Evaluate(sample_pos) * oneOverFilterIntegral;
                glm::ivec2 p = glm::ivec2{x,y} + pixel;
                if(weight<=0 || p.x < tileBounds.min.x || p.y < tileBounds.min.y || p.x >= tileBounds.max.x || p.y >= tileBounds.max.y)continue;
                FilmTilePixel& tilePixel = At(glm::ivec2{x,y} + pixel);
                tilePixel.RGB += RGB * weight;
                tilePixel.weight += weight;

            }
        }
    }

    FilmTilePixel& At(const glm::ivec2& p) {
        std::size_t tileWidth = tileBounds.max.x - tileBounds.min.x;
        return pixels[(p.y - tileBounds.min.y) * tileWidth + (p.x - tileBounds.min.x)];
    }

    const FilmTilePixel& At(const glm::ivec2& p) const {
        std::size_t tileWidth = tileBounds.max.x - tileBounds.min.x;
        return pixels[(p.y - tileBounds.min.y) * tileWidth + (p.x - tileBounds.min.x)];
    }

    Bounds2i PixelBounds() const {
        return pixelBounds;
    }

    Bounds2i TileBounds() const {
        return tileBounds;
    }

private:
    friend class Film;
    Bounds2i pixelBounds;
    Bounds2i tileBounds;
    std::shared_ptr<Filter> filter;
    glm::ivec2 radius;
    double oneOverFilterIntegral;
    std::vector<FilmTilePixel> pixels;
};

class Film {
public:
    Film(const glm::ivec2& resolution, const std::shared_ptr<Filter>& filter,double maxComponentValue = std::numeric_limits<double>::infinity()) : xResolution(resolution.x), yResolution(resolution.y), screen(resolution.x*resolution.y), filter(filter), maxComponent{maxComponentValue} {

    }

    FilmTile GetFilmTile(const Bounds2i& bounds) const {
        glm::ivec2 lowerLeft = glm::clamp(bounds.min - glm::ivec2(glm::ceil(filter->Radius() - glm::vec2{0.5,0.5})),{0,0},{xResolution,yResolution}); 
        glm::ivec2 upperRight = glm::clamp(bounds.max + glm::ivec2(glm::ceil(filter->Radius() - glm::vec2{0.5,0.5})),{0,0},{xResolution,yResolution}); 
        Bounds2i tileBounds = {lowerLeft,upperRight};
        return FilmTile(bounds,tileBounds,filter);
    }

    void Merge(const FilmTile& tile) {
        std::size_t tileWidth = tile.tileBounds.max.x - tile.tileBounds.min.x;
        for(int y = tile.tileBounds.min.y;y<tile.tileBounds.max.y;y++){
            for(int x = tile.tileBounds.min.x;x<tile.tileBounds.max.x;x++){
                FilmTilePixel pixel = tile.At({x,y});
                screen[y*xResolution + x].Add(pixel);
            }
        }
    }

    void WriteImage(const std::string& filename) const {
        std::ofstream out(filename + ".ppm", std::ios::binary );
        out << "P6\n" << xResolution << " " << yResolution << "\n255\n";
        for(int i = yResolution-1;i>=0;i--){
            for(int j = 0;j<xResolution;j++){
                glm::dvec3 color = glm::dvec3{screen[i*xResolution + j].r.load(),screen[i*xResolution + j].g.load(),screen[i*xResolution + j].b.load()} / screen[i*xResolution + j].weight.load();
                color = reinhard_jodie(color);
                double r = linear_to_sRGB(color.r);
                double g = linear_to_sRGB(color.g);
                double b = linear_to_sRGB(color.b);
                out << (char)(255.999*std::max(0.0,std::min(1.0,r))) << (char)(255.999*std::max(0.,std::min(1.,g))) << (char)(255.999*std::max(0.,std::min(1.,b)));
            }
        }
        out.close();
    }

    void ReadImage(const std::string& filename) {
        //implement
    }

    glm::ivec2 Resolution() const {
        return {xResolution,yResolution};
    }
private:  
    struct AtomicPixel{
        AtomicPixel() {
            this->r.store(0);
            this->g.store(0);
            this->b.store(0);
            this->weight.store(0);
        }

        AtomicPixel(double r, double g, double b, double weight) : r(r), g(g), b(b), weight(weight){

        }

        std::atomic<double> r{0.0};
        std::atomic<double> g{0.0};
        std::atomic<double> b{0.0};
        std::atomic<double> weight{0};

        void Add(const glm::dvec3& RGB , double weight){
            r.fetch_add(RGB.r,std::memory_order_relaxed);
            g.fetch_add(RGB.g,std::memory_order_relaxed);
            b.fetch_add(RGB.b,std::memory_order_relaxed);
            this->weight.fetch_add(weight,std::memory_order_relaxed);
        }

        void Add(const FilmTilePixel& pixel){
            Add(pixel.RGB,pixel.weight);
        }

        AtomicPixel(const AtomicPixel& other) : r(other.r.load()), g(other.g.load()), b(other.b.load()), weight(other.weight.load()) {}

        AtomicPixel& operator=(const AtomicPixel& other) {
            if(other != *this){
                r = other.r.load();
                g = other.g.load();
                b = other.b.load();
                weight = other.weight.load();
            }
            return *this;
        }

        bool operator==(const AtomicPixel& other) const {
            return r.load() == other.r.load() && g.load() == other.g.load() && b.load() == other.b.load() && weight.load() == other.weight.load();
        }

        bool operator!=(const AtomicPixel& other) const {
            return !(other == *this);
        }
    };
    uint16_t xResolution;
    uint16_t yResolution;
    std::vector<AtomicPixel> screen; 
    std::shared_ptr<Filter> filter;
    double maxComponent;
};