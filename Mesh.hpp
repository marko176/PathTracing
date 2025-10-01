#pragma once
#include <vector>
#include <memory>
#include <glm/glm.hpp>
#include <cstring>
#include "Shape.hpp"
#include "Material.hpp"

//resource manager should hold shared_ptr<mesh>
//we have set<shared_ptr>
//when we create new mesh in model we check if it exists -> !

class Mesh{
public:
    Mesh(const std::vector<uint32_t>& indices,const std::vector<glm::vec3>& vertices,const std::vector<glm::vec3>& tangents,const std::vector<glm::vec3>& bitangents,const std::vector<glm::vec3>& normals, const std::vector<glm::vec2>& texCoords, const std::shared_ptr<Material>& mat) ;
    
    uint32_t GetTriangleCount() const {
        return triangle_count;
    }

    uint32_t GetVertexCount() const {
        return vertex_count;
    }

    std::vector<uint32_t> GetIndices() const {
        return indices;
    }

    std::vector<glm::vec3> GetVertices() const {
        return vertices;
    }

    std::vector<glm::vec3> GetTangents() const {
        return tangents;
    }

    std::vector<glm::vec3> GetBitangents() const {
        return bitangents;
    }

    std::vector<glm::vec3> GetNormals() const {
        return normals;
    }

    std::vector<glm::vec2> GetTexCoords() const {
        return texCoords;
    }

    std::shared_ptr<Material> GetMaterial() const {
        return material;
    }

    const std::shared_ptr<std::vector<TriangleShape>>& GetControlPtr() const {
        return shapes;
    }
    TriangleShape* GetShape(std::size_t index) const {
        return &(*shapes)[index];
    }

    bool operator==(const Mesh& other) const = default;
    bool operator!=(const Mesh& other) const {
        return !(*this == other);
    }
private:
    uint32_t triangle_count;
    uint32_t vertex_count;
    std::vector<uint32_t> indices;//store indices in glm::vec3i ?? 
    std::vector<glm::vec3> vertices;
    std::vector<glm::vec3> tangents;
    std::vector<glm::vec3> bitangents;
    std::vector<glm::vec3> normals;
    std::vector<glm::vec2> texCoords;
    std::shared_ptr<Material> material;
    std::shared_ptr<std::vector<TriangleShape>> shapes;
};

