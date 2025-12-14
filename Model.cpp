#include "Model.hpp"
#include <assimp/scene.h>
#include <assimp/Importer.hpp>
#include <assimp/Exporter.hpp>
#include <assimp/postprocess.h>
#include <assimp/material.h>
#include <assimp/GltfMaterial.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>
#include "Material.hpp"
#include "Medium.hpp"
#include "ResourceManager.hpp"


Model::Model(const std::string& path){
    if(!load_model(path)){
        std::cerr << "Failed to load model: " + path << std::endl;
    }
}

bool Model::load_model(const std::string& path){
    Assimp::Importer importer;

    const aiScene* scene = nullptr;
    model_path = GetModelDirectory(path);
    std::string format = GetFormat(path);//optional

    if(format == "assbin"){
        scene = importer.ReadFile(path, 0);
        if(!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode){
            std::cout << "ERROR::ASSIMP::" << importer.GetErrorString() << "\n";
            importer.FreeScene();
            return false;
        }
    } else if(!format.empty()){
        //tangents are still bad in some models? -> need to use mikkTspace?
        //quads cause problems? -> look for models with quads (we should use quad shapes if possible)
        scene = importer.ReadFile(path, 
            aiProcess_Triangulate |
            //aiProcess_RemoveComponent    |
            aiProcess_JoinIdenticalVertices |
            aiProcess_GenUVCoords |
            aiProcess_SortByPType |
            aiProcess_RemoveRedundantMaterials |
            aiProcess_FindInvalidData |
            aiProcess_FlipUVs |
            aiProcess_CalcTangentSpace |
            aiProcess_GenSmoothNormals |
            aiProcess_ImproveCacheLocality |
            aiProcess_OptimizeMeshes |
            aiProcess_SplitLargeMeshes |
            aiProcess_PreTransformVertices //WILL NEED TO REMOVE THIS AND TO TRANSFORMATIONS OUR SELF
            //aiProcess_ValidateDataStructure
        );
        if(!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode){
            std::cout << "ERROR::ASSIMP::" << importer.GetErrorString() << "\n";
            importer.FreeScene();
            return false;
        }

        std::string outputPath = GetModelPath(path);
        Assimp::Exporter exporter;
        exporter.Export(scene, "assbin", outputPath + ".assbin");//model_path + temp.assbin
    }

    process_node(scene->mRootNode, scene);
    importer.FreeScene();
    std::cout << "Model " << path.substr(path.find_last_of('/') + 1) << " Created" << std::endl;
    return true;
}
void Model::process_node(aiNode* node, const aiScene* scene){

    for(unsigned int i = 0;i < node->mNumMeshes;i++){
        aiMesh* mesh = scene->mMeshes[node->mMeshes[i]];
        meshes.push_back(process_mesh(mesh, scene));
    }

    for(unsigned int i = 0;i < node->mNumChildren;i++){
        process_node(node->mChildren[i], scene);
    }

}

std::vector<std::shared_ptr<Texture>> Model::GetTextures(aiMaterial* material) const{
    std::vector<std::shared_ptr<Texture>> textures(AI_TEXTURE_TYPE_MAX, nullptr);
    std::vector<glm::vec3> scales(AI_TEXTURE_TYPE_MAX, glm::vec3(1));
    float intensity = 1;
    material->Get(AI_MATKEY_EMISSIVE_INTENSITY, intensity);

    std::cout << "I:" << intensity << " \n";

    scales[aiTextureType_EMISSIVE] = glm::vec3(intensity);

    aiColor3D emissiveColor(0, 0, 0);
    material->Get(AI_MATKEY_COLOR_EMISSIVE, emissiveColor);
    std::cout << "emissiveColor:" << emissiveColor.r << " " << emissiveColor.g << " " << emissiveColor.b << std::endl;
    aiColor4D baseColor(1, 1, 1, 1);
    material->Get(AI_MATKEY_BASE_COLOR, baseColor);
    std::cout << "BASE COLOR:" << baseColor.r << " " << baseColor.g << " " << baseColor.b << " " << baseColor.a << std::endl;

    scales[aiTextureType_BASE_COLOR] = glm::vec3(baseColor.r, baseColor.g, baseColor.b);
    scales[aiTextureType_DIFFUSE] = glm::vec3(baseColor.r, baseColor.g, baseColor.b);
    if(glm::vec3(emissiveColor.r, emissiveColor.g, emissiveColor.b) != glm::vec3(0) && scales[aiTextureType_EMISSIVE] != glm::vec3(0))textures[aiTextureType_EMISSIVE] = std::make_shared<SolidColor>(glm::vec3(emissiveColor.r, emissiveColor.g, emissiveColor.b), scales[aiTextureType_EMISSIVE]);
    //create image here and do convolution with baseColor
    //create texture(image)

    for(int i = 0;i < AI_TEXTURE_TYPE_MAX;i++){
        aiString texturePath;
        material->GetTexture(aiTextureType(i), 0, &texturePath);
        std::string path = model_path + texturePath.C_Str();
        std::replace(path.begin(), path.end(), '\\', '/');

        if(path != model_path){
            std::cout << i << " " << path << "\n";
            std::shared_ptr<Texture> tex = ResourceManager::get_instance().GetImageTexture(path, i == aiTextureType_DIFFUSE || i == aiTextureType_BASE_COLOR || i == aiTextureType_EMISSIVE, scales[i]);//gamma corection when i == i == aiTextureType_DIFFUSE
            textures[i] = tex;
        }
    }
    return textures;
}

std::shared_ptr<Material> Model::SetupMaterial(const std::vector<std::shared_ptr<Texture>>& textures, aiMaterial* material) const{
    std::shared_ptr<Material> mat = nullptr;
    float ior = 1.5;
    material->Get(AI_MATKEY_REFRACTI, ior);
    std::cout << "ior: " << ior << std::endl;
    float opacity = 1.0f;
    material->Get(AI_MATKEY_OPACITY, opacity);
    std::cout << "opacity: " << opacity << std::endl;
    float transmission = 0.0f;
    bool hasTransmission = material->Get(AI_MATKEY_TRANSMISSION_FACTOR, transmission) == AI_SUCCESS;
    std::cout << "transmission: " << transmission << std::endl;
    float thickeness = 0;
    bool IsThinDielectric = material->Get(AI_MATKEY_VOLUME_THICKNESS_FACTOR, thickeness) == AI_SUCCESS && thickeness == 0;

    std::cout << "thickeness: " << thickeness << std::endl;
    aiColor4D baseColor(1, 1, 1, 1);
    material->Get(AI_MATKEY_BASE_COLOR, baseColor);

    aiColor3D transmissionColor(1, 1, 1);
    material->Get(AI_MATKEY_COLOR_TRANSPARENT, transmissionColor);
    std::cout << "transmissionColor COLOR:" << transmissionColor.r << " " << transmissionColor.g << " " << transmissionColor.b << " " << std::endl;

    float metallicFactor = 1;
    bool hasMetallicFactor = material->Get(AI_MATKEY_METALLIC_FACTOR, metallicFactor) == AI_SUCCESS;
    std::cout << "metallicFactor:" << metallicFactor << " " << hasMetallicFactor << std::endl;

    AlphaTester alphaTester {};
    aiString alphaMode;

    if(material->Get(AI_MATKEY_GLTF_ALPHAMODE, alphaMode) == AI_SUCCESS){
        std::string tmp = alphaMode.C_Str();
        std::cout << tmp << std::endl;
        if(tmp == "OPAQUE"){
            alphaTester.mode = AlphaMode::Opaque;
        } else if(tmp == "BLEND"){
            alphaTester.mode = AlphaMode::Blend;
        } else if(tmp == "MASK"){
            alphaTester.mode = AlphaMode::Mask;
            material->Get(AI_MATKEY_GLTF_ALPHACUTOFF, alphaTester.cutoff);
        }
    }

    //we should just add the alpha texture values to the rgb tex values

    std::shared_ptr<Texture> albedo = textures[aiTextureType_BASE_COLOR] != nullptr ? textures[aiTextureType_BASE_COLOR] : textures[aiTextureType_DIFFUSE];
    if(albedo == nullptr)albedo = std::make_shared<SolidColor>(glm::vec3 { baseColor.r,baseColor.g,baseColor.b });
    std::shared_ptr<Texture> normal = textures[aiTextureType_NORMALS] != nullptr ? textures[aiTextureType_NORMALS] : textures[aiTextureType_HEIGHT];
    std::shared_ptr<Texture> roughness = textures[aiTextureType_DIFFUSE_ROUGHNESS] != nullptr ? textures[aiTextureType_DIFFUSE_ROUGHNESS] : textures[aiTextureType_SHININESS];
    std::shared_ptr<Texture> metallic = textures[aiTextureType_METALNESS] != nullptr ? textures[aiTextureType_METALNESS] : textures[aiTextureType_AMBIENT];
    std::shared_ptr<Texture> alpha = textures[aiTextureType_OPACITY];//fix this to support base alpha

    //alpha texture max 0b1000
    //if we have 4 channels ? stbi_info to get channel count

    if(hasTransmission){//hasIor doesnt work for obj
        //dielectric ?
        if(IsThinDielectric){
            return std::make_shared<ThinDielectric>(ior, albedo);
        } else{
            auto tmp = std::make_shared<MicrofacetDielectric>(ior, albedo, normal, roughness);
            //if(textures[aiTextureType_TRANSMISSION])
            //    tmp->SetTransmissionTexture(textures[aiTextureType_TRANSMISSION]);
            tmp->setAlphaTester(alphaTester);
            return tmp;
        }
    } else if(hasMetallicFactor){
        auto tmp = std::make_shared<MicrofacetDiffuse>(albedo, normal, roughness, metallic, textures[aiTextureType_OPACITY]);
        tmp->setAlphaTester(alphaTester);//do this with other materials
        return tmp;
    }




    if(textures[aiTextureType_DIFFUSE] || textures[aiTextureType_BASE_COLOR]){

        //convolute between base color * diffuse texture
        //aiTextureType_DIFFUSE_ROUGHNESS
        //aiTextureType_DIFFUSE_ROUGHNESS
        auto tmp = std::make_shared<MicrofacetDiffuse>(albedo, normal, roughness, metallic, textures[aiTextureType_OPACITY]);
        tmp->setAlphaTester(alphaTester);//do this with other materials
        return tmp;
    } else{
        mat = std::make_shared<MicrofacetDiffuse>(glm::vec3(.65, .05, .05));

        //should be for eg galss get index of refraction

        /*
        if(thickeness != 0){
            aiColor3D baseColor(1,1,1);
            aiColor3D specularColor(1,1,1);
            aiColor3D transparentColor(1,1,1);
            aiColor3D volumeColor(1,1,1);
            float distance = 0;
            material->Get(AI_MATKEY_COLOR_DIFFUSE,baseColor);
            material->Get(AI_MATKEY_COLOR_SPECULAR,specularColor);
            material->Get(AI_MATKEY_COLOR_TRANSPARENT,transparentColor);
            material->Get(AI_MATKEY_VOLUME_ATTENUATION_COLOR,volumeColor);
            material->Get(AI_MATKEY_VOLUME_ATTENUATION_DISTANCE,distance);
            std::cout<<baseColor.r<<" "<<baseColor.r<<" "<<baseColor.g<<"\n";
            std::cout<<specularColor.r<<" "<<specularColor.r<<" "<<specularColor.g<<"\n";
            std::cout<<transparentColor.r<<" "<<transparentColor.r<<" "<<transparentColor.g<<"\n";
            std::cout<<volumeColor.r<<" "<<volumeColor.r<<" "<<volumeColor.g<<"\n";
            std::cout<<distance<<"\n";//to use this we need mesh to shape  -> geometric primitive
            return std::make_shared<MicrofacetDielectric>(ior,0,glm::vec3(baseColor.r,baseColor.g,baseColor.b));//maybe set as color kdlum ?? not the same
        }else{
            return std::make_shared<ThinDielectric>(ior,nullptr);
        }
        */


        if(opacity < 0.99f){

            mat = std::make_shared<MicrofacetDielectric>(1.5, 0, glm::vec3(1));//maybe set as color kdlum ?? not the same
        }

        if(material->Get(AI_MATKEY_REFRACTI, ior) == AI_SUCCESS){
            aiColor3D Kd(1, 1, 1), Ks(0, 0, 0);
            material->Get(AI_MATKEY_COLOR_DIFFUSE, Kd);
            material->Get(AI_MATKEY_COLOR_SPECULAR, Ks);

            // compute average luminance
            float kdLum = (Kd.r + Kd.g + Kd.b) / 3.0f;
            float ksLum = (Ks.r + Ks.g + Ks.b) / 3.0f;

            if(opacity < 0.99f){
                mat = std::make_shared<MicrofacetDielectric>(1.5, 0, glm::vec3(1));//maybe set as color kdlum ?? not the same
            } else if(ksLum > 0.1 && (/*kdLum == ksLum ||*/ ksLum >= 0.4)){
                //this is wrong but gives good results in san miguel
                float r = linear_to_sRGB(Ks.r);
                float g = linear_to_sRGB(Ks.g);
                float b = linear_to_sRGB(Ks.b);
                mat = std::make_shared<SpecularConductor>(glm::vec3(r, g, b));
            } else if(kdLum < ksLum){

                mat = std::make_shared<MicrofacetDielectric>(1.33, 0, glm::vec3(0.98, 1, 1));//maybe set as color kdlum ?? not the same
            } else if(ksLum > 0.1){
                mat = std::make_shared<MicrofacetDielectric>(1.5, 0, glm::vec3(1, 1, 1));//maybe set as color kdlum ?? not the same
            } else if(kdLum > 0.1 && ksLum < 0.03){
                //remove ksLum< 0.03
                //std::cout<<"Dielectric:" <<Kd.r << " " <<Kd.g << " "<<Kd.b << " "<< " "<<Ks.r << " " <<Ks.g << " "<<Ks.b<<"\n";
                float r = Ks.r + Kd.r;
                float g = Ks.g + Kd.g;
                float b = Ks.b + Kd.b;
                mat = std::make_shared<MicrofacetDielectric>(1.33, 0, glm::vec3(r, g, b));

            } else{
                float r = Kd.r;
                float g = Kd.g;
                float b = Kd.b;
                auto color = std::make_shared<SolidColor>(glm::vec3(r, g, b));
                mat = std::make_shared<MicrofacetDiffuse>(color);

                //std::cout<<Kd.r << " " <<Kd.g << " "<<Kd.b << " "<< " "<<ksLum<<"\n";
            }
        }

    }
    return mat;
}




//we should pass function to process mesh!
//std::shared_ptr<Mesh> foo(aiMesh* mesh, const aiScene* scene){}
std::shared_ptr<Mesh> Model::process_mesh(aiMesh* mesh, const aiScene* scene){

    unsigned int n = mesh->mNumVertices;

    std::vector<glm::vec3> vertices;
    std::vector<glm::vec3> tangents;
    std::vector<glm::vec3> normals;
    std::vector<glm::vec2> texCoords;
    std::vector<uint32_t> indices;
    vertices.reserve(n);
    tangents.reserve(n);
    normals.reserve(n);
    texCoords.reserve(n * 3 / 2);
    indices.reserve(n * 3 / 2);
    for(unsigned int i = 0;i < n;i++){
        if(mesh->HasPositions()){
            glm::vec3 pos = { mesh->mVertices[i].x,mesh->mVertices[i].y,mesh->mVertices[i].z };
            vertices.push_back(pos);
        }
        if(mesh->HasNormals()){
            glm::vec3 pos = { mesh->mNormals[i].x,mesh->mNormals[i].y,mesh->mNormals[i].z };
            normals.push_back(pos);
        }
        if(mesh->HasTangentsAndBitangents()){
            glm::vec3 tangent = { mesh->mTangents[i].x,mesh->mTangents[i].y,mesh->mTangents[i].z };
            tangents.push_back(tangent);
        }

        if(mesh->HasTextureCoords(0)){
            glm::vec2 pos = { mesh->mTextureCoords[0][i].x,mesh->mTextureCoords[0][i].y };
            texCoords.push_back(pos);
        } else{
            texCoords.push_back(glm::vec2 { 0,0 });
        }


    }

    for(unsigned int i = 0; i < mesh->mNumFaces; i++){
        aiFace face = mesh->mFaces[i];
        for(unsigned int j = 0; j < face.mNumIndices; j++)
            indices.push_back(face.mIndices[j]);

    }

    if(mesh->mMaterialIndex >= 0){
        aiMaterial* material = scene->mMaterials[mesh->mMaterialIndex];
        std::vector<std::shared_ptr<Texture>> textures = GetTextures(material);

        std::shared_ptr<Material> mat = SetupMaterial(textures, material);
        std::shared_ptr<Medium> medium = nullptr;
        float thickeness = 0;

        aiColor3D volumeColor(1, 1, 1);
        float distance = std::numeric_limits<float>::infinity();
        material->Get(AI_MATKEY_VOLUME_ATTENUATION_COLOR, volumeColor);
        material->Get(AI_MATKEY_VOLUME_ATTENUATION_DISTANCE, distance);
        std::cout << "volume color:" << volumeColor.r << " " << volumeColor.r << " " << volumeColor.g << "\n";
        std::cout << "distance:" << distance << "\n";//to use this we need mesh to shape  -> geometric primitive
        if(material->Get(AI_MATKEY_VOLUME_THICKNESS_FACTOR, thickeness) == AI_SUCCESS && thickeness > 0){
            glm::vec3 sigma_a = -glm::log(glm::vec3(volumeColor.r, volumeColor.g, volumeColor.b)) / distance;
            medium = std::make_shared<HomogeneusMedium>(sigma_a, glm::vec3 { 0 }, std::make_shared<HenyeyGreenstein>(0.0), 1);
        }
        return ResourceManager::get_instance().getMesh(indices, vertices, tangents, normals, texCoords, mat, textures[aiTextureType_EMISSIVE], medium);
    }
    std::shared_ptr<Texture> tmp = std::make_shared<SolidColor>(glm::vec3(.65, .05, .05));
    std::shared_ptr<Material> mat = std::make_shared<MicrofacetDiffuse>(tmp);
    return ResourceManager::get_instance().getMesh(indices, vertices, tangents, normals, texCoords, mat, nullptr, nullptr);
}

std::vector<std::shared_ptr<Mesh>> Model::GetMeshes() const{
    return meshes;
}