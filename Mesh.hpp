#pragma once
#include <vector>
#include <memory>
#include <glm/glm.hpp>
#include <cstring>
#include "Shape.hpp"
#include "Material.hpp"



class Mesh{
public:
    Mesh(const std::vector<uint32_t>& indices,const std::vector<glm::vec3>& vertices,const std::vector<glm::vec3>& tangents,const std::vector<glm::vec3>& bitangents,const std::vector<glm::vec3>& normals, const std::vector<glm::vec2>& texCoords, const std::shared_ptr<Material>& mat) ;
    uint32_t triangle_count;
    uint32_t vertex_count;
    std::vector<uint32_t> indices;//store indices in glm::vec3i ?? 
    std::vector<glm::vec3> vertices;
    std::vector<glm::vec3> tangents;
    std::vector<glm::vec3> bitangents;
    std::vector<glm::vec3> normals;
    std::vector<glm::vec2> texCoords;
    std::shared_ptr<Material> material;
    
    const std::shared_ptr<std::vector<TriangleShape>>& getControlPtr() const {
        return shapes;
    }
    TriangleShape* getShape(std::size_t index) const {
        return &(*shapes)[index];
    }
private:
    std::shared_ptr<std::vector<TriangleShape>> shapes;
};

