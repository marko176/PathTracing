#include "Mesh.hpp"
#include <memory>
#include <iostream>
Mesh::Mesh(const std::vector<uint32_t>& indices,const std::vector<glm::vec3>& vertices,const std::vector<glm::vec3>& tangents,const std::vector<glm::vec3>& bitangents,const std::vector<glm::vec3>& normals, const std::vector<glm::vec2>& texCoords, const std::shared_ptr<Material>& mat) : triangle_count(indices.size()/3), vertex_count(vertices.size()), indices(indices), vertices(vertices), tangents(tangents),bitangents(bitangents),normals(normals), texCoords(texCoords), material(mat) {
        
    shapes.reserve(triangle_count);
    for(int j = 0;j<triangle_count;j++){
        shapes.emplace_back(0,j);
    }
}