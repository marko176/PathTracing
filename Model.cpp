#include "Model.hpp"
#include <assimp/scene.h>
#include <assimp/Importer.hpp>
#include <assimp/Exporter.hpp>
#include <assimp/postprocess.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>
#include "Material.hpp"
#include "Medium.hpp"
#include "ResourceManager.hpp"
#include <memory_resource>

Model::Model(const std::string& path){
    if(!load_model(path)){
        std::cerr << "Failed to load model: " + path<<std::endl;
    }else{
        std::vector<GeometricPrimitive> primitives;
        primitives.reserve(1'000'000);
        //auto mat = std::make_shared<dielectric>(1.5,glm::vec3(1));
        //std::shared_ptr<Medium> med = std::make_shared<HomogeneusMedium>(glm::vec3(0.0,0,0),glm::vec3(0.01,0.9,0.9),50.0f);
        for(const std::shared_ptr<Mesh>& m : meshes){
            uint32_t n = m->triangle_count;
            for(uint32_t j = 0;j<n;j++){
                std::shared_ptr<Shape> shape = std::shared_ptr<Shape>(m->GetControlPtr(),m->GetShape(j));
                std::shared_ptr<AreaLight> area = m->GetEmissiveTexture() != nullptr ? std::make_shared<AreaLight>(shape,m->GetEmissiveTexture()) : nullptr;
                if(area){
                    area->PreProcess({});
                    if(area->Power()==0)area = nullptr;
                }
                primitives.emplace_back(shape,m->material,area,nullptr);
            }
        }
        model_bvh = BLAS(std::move(primitives));
    }
}

Model::Model(const std::string& path,const std::shared_ptr<Material>& material, const std::shared_ptr<Medium>& medium){
    if(!load_model(path)){
        std::cerr << "Failed to load model: " + path<<std::endl;
    }else{
        std::vector<GeometricPrimitive> primitives;
        primitives.reserve(1'000'000);

        for(const std::shared_ptr<Mesh>& m : meshes){
            uint32_t n = m->triangle_count;
            for(uint32_t j = 0;j<n;j++){
                std::shared_ptr<Shape> shape = std::shared_ptr<Shape>(m->GetControlPtr(),m->GetShape(j));
                std::shared_ptr<AreaLight> area = m->GetEmissiveTexture() != nullptr ? std::make_shared<AreaLight>(shape,m->GetEmissiveTexture()) : nullptr;
                if(area){
                    area->PreProcess({});
                    if(area->Power()==0)area = nullptr;//this keeps black parts of emmisive texture out // should be 0 maybe?
                }
                primitives.emplace_back(shape,material,area,medium);
            }
        }
        model_bvh = BLAS(std::move(primitives));
    }
}



bool Model::load_model(const std::string& path) {
    Assimp::Importer importer;

    const aiScene* scene = nullptr;
    model_path = GetModelDirectory(path);
    std::size_t index = path.find_last_of('.');
    Format format = GetFormat(path);

    if(format == Format::ASSBIN){
        scene = importer.ReadFile(path,0);
        if(!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode){
            std::cout << "ERROR::ASSIMP::" << importer.GetErrorString() << "\n";
            importer.FreeScene();
            return false;
        }
    }else if(format != Format::NONE){
        //need to fix degenerate triangles
        scene = importer.ReadFile(path,  aiProcess_Triangulate |
                                                    //aiProcess_RemoveComponent    |
		                                            aiProcess_JoinIdenticalVertices | //not needed?
		                                            aiProcess_GenUVCoords |
		                                            aiProcess_SortByPType |
		                                            aiProcess_RemoveRedundantMaterials |
		                                            aiProcess_FindInvalidData |
		                                            aiProcess_FlipUVs |
		                                            aiProcess_CalcTangentSpace |
		                                            aiProcess_GenSmoothNormals |
		                                            aiProcess_ImproveCacheLocality |
		                                            aiProcess_OptimizeMeshes |
		                                            aiProcess_SplitLargeMeshes 
                                                    //aiProcess_ValidateDataStructure
                                                    );
        if(!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode){
            std::cout << "ERROR::ASSIMP::" << importer.GetErrorString() << "\n";
            importer.FreeScene();
            return false;
        }

        std::string outputPath = GetModelPath(path);
        Assimp::Exporter exporter;
        exporter.Export(scene,"assbin",outputPath + ".assbin");//model_path + temp.assbin
    }

    process_node(scene->mRootNode,scene);
    importer.FreeScene();
    std::cout<<"Model "<< path.substr(path.find_last_of('/') + 1)<<" Created"<<std::endl;
    return true;
}
void Model::process_node(aiNode* node, const aiScene* scene) {

    for(unsigned int i = 0;i< node->mNumMeshes;i++){
        aiMesh* mesh = scene->mMeshes[node->mMeshes[i]];
        meshes.push_back(process_mesh(mesh,scene));
    }

    for(unsigned int i = 0;i<node->mNumChildren;i++){
        process_node(node->mChildren[i],scene);
    }

}

std::vector<std::shared_ptr<Texture>> Model::GetTextures(aiMaterial* material) const {
    std::vector<std::shared_ptr<Texture>> textures(AI_TEXTURE_TYPE_MAX);
    //textures[aiTextureType_SHININESS] = std::make_shared<SolidColor>(glm::vec3(1));//use resource manager
    //textures[aiTextureType_AMBIENT] = std::make_shared<SolidColor>(glm::vec3(0));//use resource manager
    //if nullptr then do it in material 
    for(int i = 0;i<AI_TEXTURE_TYPE_MAX;i++){
        aiString texturePath;
        material->GetTexture(aiTextureType(i),0,&texturePath);
        std::string path = model_path + texturePath.C_Str();
        std::replace(path.begin(),path.end(),'\\','/');
        std::shared_ptr<Texture> tex = nullptr;
        if(path != model_path){
            std::cout<<i<<" "<<path<<"\n";
            
            tex = ResourceManager::get_instance().GetImageTexture(path,i == aiTextureType_DIFFUSE || i == aiTextureType_BASE_COLOR);//gamma corection when i == i == aiTextureType_DIFFUSE
        }
        textures[i]=tex;
    }
    
    return textures;
}

std::shared_ptr<Material> Model::SetupOBJMaterial(const std::vector<std::shared_ptr<Texture>>& textures,aiMaterial* material) const {
    std::shared_ptr<Material> mat = nullptr;

    if(textures[aiTextureType_DIFFUSE]){
        //aiTextureType_DIFFUSE_ROUGHNESS
        //aiTextureType_DIFFUSE_ROUGHNESS
        std::shared_ptr<Texture> albedo = textures[aiTextureType_BASE_COLOR] != nullptr ? textures[aiTextureType_BASE_COLOR] : textures[aiTextureType_DIFFUSE];
        std::shared_ptr<Texture> normal = textures[aiTextureType_NORMALS] != nullptr ? textures[aiTextureType_NORMALS] : textures[aiTextureType_HEIGHT];
        std::shared_ptr<Texture> roughness = textures[aiTextureType_DIFFUSE_ROUGHNESS] != nullptr ? textures[aiTextureType_DIFFUSE_ROUGHNESS] : textures[aiTextureType_SHININESS];
        std::shared_ptr<Texture> metallic = textures[aiTextureType_METALNESS] != nullptr ? textures[aiTextureType_METALNESS] : textures[aiTextureType_AMBIENT];
        mat = std::make_shared<lambertian>(albedo,normal,roughness,metallic,textures[aiTextureType_OPACITY]);
    }else{
        
        mat = std::make_shared<lambertian>(glm::vec3(.65, .05, .05));
        //should be for eg galss get index of refraction
        float ior = 1;
        if(material->Get(AI_MATKEY_REFRACTI,ior) == AI_SUCCESS){
            aiColor3D Kd(1,1,1), Ks(0,0,0);
            material->Get(AI_MATKEY_COLOR_DIFFUSE,  Kd);
            material->Get(AI_MATKEY_COLOR_SPECULAR, Ks);

            // compute average luminance
            float kdLum = (Kd.r + Kd.g + Kd.b) / 3.0f;
            float ksLum = (Ks.r + Ks.g + Ks.b) / 3.0f;
            
            float opacity = 1.0f;
            material->Get(AI_MATKEY_OPACITY, opacity);
            if(opacity < 0.99f){
                mat = std::make_shared<MicrofacetDielectric>(1.5,0,glm::vec3(1));//maybe set as color kdlum ?? not the same
            }else if(ksLum > 0.1 && (/*kdLum == ksLum ||*/ ksLum >= 0.4)){
                //this is wrong but gives good results in san miguel
                float r = linear_to_sRGB(Ks.r);
                float g = linear_to_sRGB(Ks.g);
                float b = linear_to_sRGB(Ks.b);
                mat = std::make_shared<SpecularConductor>(glm::vec3(r,g,b));
            }else if(kdLum < ksLum){
                
                mat = std::make_shared<MicrofacetDielectric>(1.33,0,glm::vec3(0.98,1,1));//maybe set as color kdlum ?? not the same
            }else if(ksLum > 0.1){
                mat = std::make_shared<MicrofacetDielectric>(1.5,0,glm::vec3(1,1,1));//maybe set as color kdlum ?? not the same
            }else if(kdLum > 0.1 && ksLum < 0.03){
                //remove ksLum< 0.03
                //std::cout<<"Dielectric:" <<Kd.r << " " <<Kd.g << " "<<Kd.b << " "<< " "<<Ks.r << " " <<Ks.g << " "<<Ks.b<<"\n";
                float r = Ks.r + Kd.r;
                float g = Ks.g + Kd.g;
                float b = Ks.b + Kd.b;
                mat = std::make_shared<MicrofacetDielectric>(1.33,0,glm::vec3(r,g,b));
                
            }else{
                float r = Kd.r;
                float g = Kd.g;
                float b = Kd.b;
                auto color = std::make_shared<SolidColor>(glm::vec3(r,g,b));
                mat = std::make_shared<lambertian>(color);
    
                //std::cout<<Kd.r << " " <<Kd.g << " "<<Kd.b << " "<< " "<<ksLum<<"\n";
            }
        }
        
    }
    return mat;
}




//we should pass function to process mesh!
//std::shared_ptr<Mesh> foo(aiMesh* mesh, const aiScene* scene){}
std::shared_ptr<Mesh> Model::process_mesh(aiMesh* mesh, const aiScene* scene) {

    unsigned int n = mesh->mNumVertices;

    std::vector<glm::vec3> vertices;
    std::vector<glm::vec3> tangents;
    std::vector<glm::vec3> bitangents;
    std::vector<glm::vec3> normals;
    std::vector<glm::vec2> texCoords;
    std::vector<uint32_t> indices;
    vertices.reserve(n);
    tangents.reserve(n);
    bitangents.reserve(n);
    normals.reserve(n);
    texCoords.reserve(n*3/2);
    indices.reserve(n*3/2);
    for(unsigned int i = 0;i<n;i++){
    
        if(mesh->HasPositions()){
            glm::vec3 pos = {mesh->mVertices[i].x,mesh->mVertices[i].y,mesh->mVertices[i].z};
            vertices.push_back(pos);
        }
        if(mesh->HasNormals()){
            glm::vec3 pos = {mesh->mNormals[i].x,mesh->mNormals[i].y,mesh->mNormals[i].z};
            normals.push_back(pos);
        }
        if(mesh->HasTangentsAndBitangents()){
            glm::vec3 tangent = {mesh->mTangents[i].x,mesh->mTangents[i].y,mesh->mTangents[i].z};
            tangents.push_back(tangent);
           
            glm::vec3 bitangent = {mesh->mBitangents[i].x,mesh->mBitangents[i].y,mesh->mBitangents[i].z};
            bitangents.push_back(bitangent);
        }

        if(mesh->HasTextureCoords(0)){
            glm::vec2 pos = {mesh->mTextureCoords[0][i].x,mesh->mTextureCoords[0][i].y};
            texCoords.push_back(pos);
        }else{
            texCoords.push_back(glm::vec2{0,0});
        }

    
    }

    for(unsigned int i = 0; i < mesh->mNumFaces; i++)
    {
        aiFace face = mesh->mFaces[i];
        for(unsigned int j = 0; j < face.mNumIndices; j++)
            indices.push_back(face.mIndices[j]);
        
    }  

    //if(mesh->mMaterialIndex >= 0) setupMaterial(scene->mMaterials[mesh->mMaterialIndex],obj/gltf )

    //if is emmisive we need light which holds prt to this mesh?
    //
    if(mesh->mMaterialIndex >= 0){
        aiMaterial* material = scene->mMaterials[mesh->mMaterialIndex];
        std::vector<std::shared_ptr<Texture>> textures = GetTextures(material);

        std::shared_ptr<Material> mat = SetupOBJMaterial(textures,material);
        //std::shared_ptr<Texture> tmp = std::make_shared<SolidColor>(glm::vec3(.65, .05, .05));
        //std::shared_ptr<Material> mat = std::make_shared<lambertian>(tmp);
        return ResourceManager::get_instance().getMesh(indices,vertices,tangents,bitangents,normals,texCoords,mat,textures[aiTextureType_EMISSIVE]);
        //return std::make_shared<Mesh>(indices,vertices,tangents,bitangents,normals,texCoords,mat);
    }
    std::shared_ptr<Texture> tmp = std::make_shared<SolidColor>(glm::vec3(.65, .05, .05));
    std::shared_ptr<Material> mat = std::make_shared<lambertian>(tmp);
    return ResourceManager::get_instance().getMesh(indices,vertices,tangents,bitangents,normals,texCoords,mat,nullptr);
    //return std::make_shared<Mesh>(indices,vertices,tangents,bitangents,normals,texCoords,mat);
}

auto Model::GetMeshes() const -> const std::vector<std::shared_ptr<Mesh>>& {
    return meshes;
}