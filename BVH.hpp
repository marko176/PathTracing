#pragma once
#include "Primitive.hpp"
#include <atomic>
#include <thread>
#include <execution>
inline constexpr unsigned char P(int a, int b, int c, int d){
    return (a << 6) | (b << 4) | (c << 2) | d;
}

inline constexpr unsigned char PermToIndexLUT[24] = {
    P(0,1,2,3), P(0,1,3,2), P(0,2,1,3), P(0,2,3,1),
    P(0,3,1,2), P(0,3,2,1), P(1,0,2,3), P(1,0,3,2),
    P(1,2,0,3), P(1,2,3,0), P(1,3,0,2), P(1,3,2,0),
    P(2,0,1,3), P(2,0,3,1), P(2,1,0,3), P(2,1,3,0),
    P(2,3,0,1), P(2,3,1,0), P(3,0,1,2), P(3,0,2,1),
    P(3,1,0,2), P(3,1,2,0), P(3,2,0,1), P(3,2,1,0)
};

inline constexpr unsigned char PermToIndex(unsigned int perm){
    for(int i = 0;i < 24;i++){
        if(perm == PermToIndexLUT[i])return i;
    }
    return 0;
}

//BVH 4
//

struct BVH_NODE{
    AABB bbox;
    uint32_t left = 0;
    uint32_t right = 0;// right node / first triangle / meshID  / modelID
    uint16_t count = 0;// is_leaf / triangle count / mesh count / modelCount
    uint16_t axis = 0;  //holds split axis, left is - of axis, right is + of axis
    //is ray sign is + -> left is closer (push right then left)
    bool isLeaf() const{
        return count != 0;
    }
};

struct BVH4_NODE{
    unsigned char count = 0;//if leaf this holds triangle count
    unsigned char active = 0b0000;//if it is 0 -> leaf node
    unsigned char perm = 0;//build bvh( parent node, child index) -> parent[node][child index] = perm
    uint32_t ClusterIdx = 0;//if leaf then this holds triangle index
};

struct alignas(128) BVH4_CLUSTER{
    float xmin[4] { 0 }, xmax[4] { 0 };
    float ymin[4] { 0 }, ymax[4] { 0 };
    float zmin[4] { 0 }, zmax[4] { 0 };

    BVH4_NODE children[4] {};

    void SetDataAtIndex(const AABB& bbox, int index){
        xmin[index] = bbox.min[0];
        ymin[index] = bbox.min[1];
        zmin[index] = bbox.min[2];
        xmax[index] = bbox.max[0];
        ymax[index] = bbox.max[1];
        zmax[index] = bbox.max[2];
    }
};



template <typename T>
class BVHBase : public Primitive{
public:
    virtual ~BVHBase() = default;
    BVHBase() = default;
    std::vector<std::shared_ptr<Light>> GetLights() const final{
        std::vector<std::shared_ptr<Light>> lights;
        for(const T& prim : primitives){
            std::vector<std::shared_ptr<Light>> primLights;
            if constexpr(requires { prim->GetLights(); }){
                primLights = prim->GetLights();
            } else{
                primLights = prim.GetLights();
            }
            lights.insert(lights.end(), primLights.begin(), primLights.end());
        }
        return lights;
    }

    AABB BoundingBox() const final{
        return BVHbbox;
    }

protected:
    struct PrimitiveInfo{
        PrimitiveInfo(uint32_t index, const AABB& bbox) : index { index }, bbox { bbox }, centroid(0.5f * (bbox.max + bbox.min)){}
        uint32_t index;
        AABB bbox;
        glm::vec3 centroid;
    };

    static std::vector<BVH_NODE> BuildBVHBase(std::vector<PrimitiveInfo>& primitiveInfo){
        //std::vector<BVH_NODE> baseNodes;
        //baseNodes.reserve(primitiveInfo.size() * 2 - 1);
        //BVHBase::BuildBase(0, primitiveInfo.size(), primitiveInfo, baseNodes);
        std::vector<BVH_NODE> baseNodes(primitiveInfo.size() * 2 - 1);
        std::atomic<uint32_t> idx { 0 };
        BVHBase::BuildBaseThreaded(0, primitiveInfo.size(), primitiveInfo, baseNodes, idx);
        baseNodes.resize(idx.load(std::memory_order_relaxed));
        return baseNodes;
    }

    inline bool intersectPrimitives(const Ray& ray, float& max, uint32_t primIndex, uint32_t primCount, SurfaceInteraction& interaction) const{
        bool hit = false;
        for(uint32_t i = primIndex;i < primIndex + primCount;i++){
            if constexpr(requires { primitives[i]->Intersect(ray, interaction, max); }){
                if(primitives[i]->Intersect(ray, interaction, max)){
                    hit = true;
                    max = interaction.t;
                }
            } else{
                if(primitives[i].Intersect(ray, interaction, max)){
                    hit = true;
                    max = interaction.t;
                }
            }
        }

        return hit;
    }

    inline bool intersectPrimitivesPred(const Ray& ray, float max, uint32_t primIndex, uint32_t primCount) const{
        for(uint32_t i = primIndex;i < primIndex + primCount;i++){
            if constexpr(requires { primitives[i]->IntersectPred(ray, max); }){
                if(primitives[i]->IntersectPred(ray, max))return true;
            } else{
                if(primitives[i].IntersectPred(ray, max))return true;
            }
        }
        return false;
    }
private:
    static uint32_t BuildBase(uint32_t first_triangle, uint32_t last_triangle, std::vector<PrimitiveInfo>& primitiveInfo, std::vector<BVH_NODE>& nodes){
        uint32_t index = nodes.size();
        BVH_NODE& node = nodes.emplace_back();
        uint32_t object_span = last_triangle - first_triangle;

        glm::vec3 maxAxis { -std::numeric_limits<float>::max() };
        glm::vec3 minAxis { std::numeric_limits<float>::max() };
        for(uint32_t i = first_triangle;i < last_triangle;i++){
            node.bbox.Expand(primitiveInfo[i].bbox);
            glm::vec3 triangle_centroid = primitiveInfo[i].centroid;
            maxAxis = glm::max(maxAxis, triangle_centroid);
            minAxis = glm::min(minAxis, triangle_centroid);
        }

        node.right = first_triangle;//triangle
        node.count = object_span;
        if(object_span > 2){

            unsigned char best_axis = 0;
            float bestPos = 0;
            float bestCost = std::numeric_limits<float>::max();
            uint32_t BINS = object_span >= 1024 ? 32 :
                object_span >= 64 ? 16 : 8;
            std::vector<Bin> bins;
            std::vector<float> rightArea;
            bins.reserve(BINS);
            rightArea.reserve(BINS - 1);



            for(unsigned char axis = 0;axis < 3;axis++){
                bins.assign(BINS, Bin {});



                //do this outside but vectorised?
                float max = maxAxis[axis];
                float min = minAxis[axis];
                /*
                for(int i = first_triangle;i<last_triangle;i++){
                    float triangle_centroid = primitiveInfo[i].centroid[axis];
                    max = std::max(max,triangle_centroid);
                    min = std::min(min,triangle_centroid);
                }
                */


                if(std::abs(max - min) <= std::numeric_limits<float>::epsilon())continue;
                float scale = BINS / (max - min);

                for(uint32_t i = first_triangle;i < last_triangle;i++){
                    //float triangle_centroid = (vertices[indices[i*3+0]][axis]+vertices[indices[i*3+1]][axis]+vertices[indices[i*3+2]][axis])/3.0f;
                    float triangle_centroid = primitiveInfo[i].centroid[axis];
                    int binId = std::min<int>(BINS - 1, (triangle_centroid - min) * scale);
                    ++bins[binId].triCount;
                    bins[binId].aabb.Expand(primitiveInfo[i].bbox);
                }

                AABB leftBox, rightBox;
                for(uint32_t i = BINS - 1;i >= 1;i--){
                    rightBox.Expand(bins[i].aabb);// i == 0 -> BINS-1  i == BINS-2 -> 1
                    rightArea[i - 1] = rightBox.Area();// i == 0 -> BINS - 2  i == BINS - 2 -> 0
                }
                /*for(int i = 0;i < BINS - 1;i++){
                    rightBox.Expand(bins[BINS - 1 - i].aabb);// i == 0 -> BINS-1  i == BINS-2 -> 1
                    rightArea[BINS - 2 - i] = rightBox.Area();// i == 0 -> BINS - 2  i == BINS - 2 -> 0
                }*/

                scale = (max - min) / BINS;
                uint32_t leftSum = 0;
                for(uint32_t i = 0;i < BINS - 1;i++){
                    leftSum += bins[i].triCount;
                    leftBox.Expand(bins[i].aabb);
                    float cost = leftSum * leftBox.Area() + (object_span - leftSum) * rightArea[i];

                    if(cost < bestCost){
                        best_axis = axis;
                        bestPos = min + (i + 1) * scale;
                        bestCost = cost;
                    }
                }

            }
            float parent_cost = node.bbox.Area() * object_span;
            if(bestCost >= parent_cost){
                return index;
            }

            uint32_t mid = std::partition(primitiveInfo.begin() + first_triangle, primitiveInfo.begin() + last_triangle, [&](const PrimitiveInfo& info){
                float triangle_centroid = info.centroid[best_axis];
                return triangle_centroid <= bestPos;
                }) - primitiveInfo.begin();


            if(mid == first_triangle || mid == last_triangle){
                return index;
            }



            //left < best pos
            //right > best pos
            //left  < right 
            //if ray sign is + we go left then right

            //node.left = pos+1;
            node.left = BuildBase(first_triangle, mid, primitiveInfo, nodes);
            node.right = BuildBase(mid, last_triangle, primitiveInfo, nodes);
            node.count = 0;
            node.axis = best_axis;

        }
        return index;
    }

    static uint32_t BuildBaseThreaded(uint32_t first_triangle, uint32_t last_triangle, std::vector<PrimitiveInfo>& primitiveInfo, std::vector<BVH_NODE>& nodes, std::atomic<uint32_t>& idx){
        uint32_t index = idx.fetch_add(1, std::memory_order_relaxed);
        BVH_NODE& node = nodes[index];
        uint32_t object_span = last_triangle - first_triangle;

        glm::vec3 maxAxis { -std::numeric_limits<float>::max() };
        glm::vec3 minAxis { std::numeric_limits<float>::max() };
        std::for_each(primitiveInfo.begin() + first_triangle, primitiveInfo.begin() + last_triangle, [&](const auto& prim){
            node.bbox.Expand(prim.bbox);
            glm::vec3 triangle_centroid = prim.centroid;
            maxAxis = glm::max(maxAxis, triangle_centroid);
            minAxis = glm::min(minAxis, triangle_centroid);
            });


        node.right = first_triangle;//triangle
        node.count = object_span;
        if(object_span > 2){

            unsigned char best_axis = 0;
            float bestPos = 0;
            float bestCost = std::numeric_limits<float>::infinity();
            uint32_t BINS = object_span >= 1024 ? 32 :
                object_span >= 64 ? 16 : 8;
            std::vector<Bin> bins;
            std::vector<float> rightArea;
            bins.reserve(BINS);
            rightArea.reserve(BINS - 1);

            for(unsigned char axis = 0;axis < 3;axis++){

                float max = maxAxis[axis];
                float min = minAxis[axis];

                if(std::abs(max - min) <= std::numeric_limits<float>::epsilon())continue;
                bins.assign(BINS, Bin {});
                float scale = BINS / (max - min);

                for(uint32_t i = first_triangle;i < last_triangle;i++){
                    float triangle_centroid = primitiveInfo[i].centroid[axis];
                    int binId = std::min<int>(BINS - 1, (triangle_centroid - min) * scale);
                    ++bins[binId].triCount;
                    bins[binId].aabb.Expand(primitiveInfo[i].bbox);
                }

                AABB leftBox, rightBox;
                for(uint32_t i = BINS - 1;i >= 1;i--){
                    rightBox.Expand(bins[i].aabb);
                    rightArea[i - 1] = rightBox.Area();
                }

                scale = (max - min) / BINS;
                uint32_t leftSum = 0;
                for(uint32_t i = 0;i < BINS - 1;i++){
                    leftSum += bins[i].triCount;
                    leftBox.Expand(bins[i].aabb);
                    float cost = leftSum * leftBox.Area() + (object_span - leftSum) * rightArea[i];

                    if(cost < bestCost){
                        best_axis = axis;
                        bestPos = min + (i + 1) * scale;
                        bestCost = cost;
                    }
                }

            }
            float parent_cost = node.bbox.Area() * object_span;
            //float minCost = 1.0f/2.0f + bestCost / node.bbox.Area();
            if(bestCost >= parent_cost /*|| minCost >= object_span*/){
                return index;
            }

            uint32_t mid = std::partition(primitiveInfo.begin() + first_triangle, primitiveInfo.begin() + last_triangle, [&](const PrimitiveInfo& info){
                float triangle_centroid = info.centroid[best_axis];
                return triangle_centroid <= bestPos;
                }) - primitiveInfo.begin();


            if(mid == first_triangle || mid == last_triangle){
                return index;
            }



            if(object_span > 256 * 1024){
                std::jthread left([&](){
                    node.left = BuildBaseThreaded(first_triangle, mid, primitiveInfo, nodes, idx);
                    });
                std::jthread right([&](){
                    node.right = BuildBaseThreaded(mid, last_triangle, primitiveInfo, nodes, idx);
                    });
            } else{
                node.left = BuildBaseThreaded(first_triangle, mid, primitiveInfo, nodes, idx);
                node.right = BuildBaseThreaded(mid, last_triangle, primitiveInfo, nodes, idx);
            }
            node.count = 0;
            node.axis = best_axis;

        }
        return index;
    }

protected:
    AABB BVHbbox;
    std::vector<T> primitives;
};


template <typename T>
class BVH : public BVHBase<T>{
    using typename BVHBase<T>::PrimitiveInfo;
public:
    virtual ~BVH() = default;
    BVH() = default;

    BVH(const std::vector<T>& prims){
        auto start = std::chrono::high_resolution_clock::now();
        std::cout << "Building BVH" << std::endl;

        BVHBase<T>::primitives.reserve(prims.size());

        std::vector<PrimitiveInfo> primitiveInfo;
        primitiveInfo.reserve(prims.size());

        for(std::size_t i = 0;i < prims.size();i++){
            if constexpr(requires { prims[i]->BoundingBox(); }){
                AABB bbox = prims[i]->BoundingBox();
                primitiveInfo.emplace_back(i, bbox);
            } else{
                AABB bbox = prims[i].BoundingBox();
                primitiveInfo.emplace_back(i, bbox);
            }
        }


        //nodes = BVHBase<T>::BuildBVHBase(primitiveInfo);
        std::vector<BVH_NODE> BVHnodes = BVHBase<T>::BuildBVHBase(primitiveInfo);
        nodes.reserve(BVHnodes.size());
        if(BVHnodes.empty()){
            BVHBase<T>::BVHbbox = AABB {};
        } else{
            BVHBase<T>::BVHbbox = BVHnodes[0].bbox;
        }
        BuildBVH2(BVHnodes);

        this->nodes.shrink_to_fit();
        for(const auto& info : primitiveInfo){
            BVHBase<T>::primitives.emplace_back(prims[info.index]);
        }

        auto duration = std::chrono::high_resolution_clock::now() - start;
        std::cout << "BVH built in: " << std::chrono::duration_cast<std::chrono::milliseconds>(duration) << std::endl;

        if(nodes.empty()){
            BVHBase<T>::BVHbbox = AABB {};
        } else{
            BVHBase<T>::BVHbbox = nodes[0].bbox;
        }
    }

    struct alignas(32) BVH2_NODE{
        AABB bbox;
        uint32_t right = 0;// right node / first triangle / meshID  / modelID
        uint16_t count = 0;// is_leaf / tirangle count / mesh count / modelCount
        uint16_t axis = 0;  //holds split axis, left is - of axis, right is + of axis
        //is ray sign is + -> left is closer (push right then left)
        bool isLeaf() const{
            return count != 0;
        }
    };

    uint32_t BuildBVH2(const std::vector<BVH_NODE>& BVHnodes, uint32_t nodeIdx = 0){
        const BVH_NODE& n = BVHnodes[nodeIdx];
        uint32_t index = nodes.size();
        BVH2_NODE& node = nodes.emplace_back();
        node.bbox = n.bbox;
        node.right = n.right;
        node.count = n.count;
        node.axis = n.axis;
        if(n.isLeaf())return index;
        BuildBVH2(BVHnodes, n.left);
        node.right = BuildBVH2(BVHnodes, n.right);
        return index;
    }


    bool IntersectPred(const Ray& ray, float max = std::numeric_limits<float>::infinity()) const final{
        if(!BVHBase<T>::BVHbbox.Hit(ray, max))return false;
        uint32_t stack[32];
        int i = 0;
        stack[i++] = 0;

        while(i){
            uint32_t index = stack[--i];
            BVH2_NODE node = nodes[index];

            // switch to intercesion test

            if(node.count == 0){
                uint32_t child1 = index + 1;
                uint32_t child2 = node.right;
                if(nodes[child1].bbox.Hit(ray, max))stack[i++] = child1;
                if(nodes[child2].bbox.Hit(ray, max))stack[i++] = child2;
            } else if(BVHBase<T>::intersectPrimitivesPred(ray, max, node.right, node.count)){
                return true;
            }

        }

        return false;
    }


    bool Intersect(const Ray& ray, SurfaceInteraction& interaction, float max = std::numeric_limits<float>::infinity()) const final{
        if(!BVHBase<T>::BVHbbox.Hit(ray, max))return false;
        uint32_t stack[32];
        int i = 0;
        stack[i++] = 0;

        bool hit_anything = false;
        bool signs[3] = { ray.dir.x > 0,ray.dir.y > 0,ray.dir.z > 0 };
        while(i){
            uint32_t index = stack[--i];
            BVH2_NODE node = nodes[index];

            if(node.count == 0){
                uint32_t child1 = index + 1;
                uint32_t child2 = node.right;
                float dist_1 = nodes[child1].bbox.HitDistance(ray, max);
                float dist_2 = nodes[child2].bbox.HitDistance(ray, max);
                uint16_t axis = node.axis;
                if(signs[axis]){
                    //if going from - -> +
                    //need to check child1 first
                    if(dist_2 != std::numeric_limits<float>::infinity())stack[i++] = child2;
                    if(dist_1 != std::numeric_limits<float>::infinity())stack[i++] = child1;
                } else{
                    if(dist_1 != std::numeric_limits<float>::infinity())stack[i++] = child1;
                    if(dist_2 != std::numeric_limits<float>::infinity())stack[i++] = child2;
                }
                /*
                if(dist_1 > dist_2){
                    std::swap(dist_1,dist_2);
                    std::swap(child1,child2);
                }
                if(dist_2 != std::numeric_limits<float>::infinity()){
                    stack[i++]=child2;
                    stack[i++]=child1;
                }else if(dist_1 != std::numeric_limits<float>::infinity()){
                    stack[i++]=child1;
                }
                */
            } else if(BVHBase<T>::intersectPrimitives(ray, max, node.right, node.count, interaction)){
                hit_anything = true;
            }

        }
        return hit_anything;
    }

private:
    std::vector<BVH2_NODE> nodes;
};

using BLAS = BVH<GeometricPrimitive>;
using TLAS = BVH<std::shared_ptr<Primitive>>;


#if defined(__SSE__) || defined(_M_AMD64) || defined(_M_X64)
template <typename T>
class BVH4 : public BVHBase<T>{
    using typename BVHBase<T>::PrimitiveInfo;
    static constexpr std::array<std::array<unsigned char, 135>, 8> LUT = [](){
        std::array<std::array<unsigned char, 135>, 8> tempLUT;
        for(unsigned int rs = 0;rs < 8;rs++){
            const unsigned int signs[3] = { rs & 1, (rs >> 1) & 1, (rs >> 2) };

            //x<0 == true -> 1
            //symetric topology 0
            for(unsigned int s0 = 0;s0 < 3;s0++){// x,y,z
                for(unsigned int s1 = 0;s1 < 3;s1++){
                    for(unsigned int s2 = 0;s2 < 3;s2++){
                        unsigned int left = 0;
                        unsigned int right = 0;
                        unsigned int perm = 0;
                        if(signs[s1]){//x < 0
                            left = 0b0100;// 1 -> 0
                        } else{
                            left = 0b0001;// 0 -> 1
                        }

                        if(signs[s2]){
                            right = 0b1110;// 3 -> 2
                        } else{
                            right = 0b1011;// 2 -> 3
                        }

                        if(signs[s0]){
                            perm = (right << 4) + left;
                        } else{
                            perm = (left << 4) + right;
                        }
                        tempLUT[rs][s0 + s1 * 3 + s2 * 9 + 0 * 27] = PermToIndex(perm);
                    }
                }
            }
            //unsymetric topology 1  -> (1 + (1 + 2))
            for(unsigned int s0 = 0;s0 < 3;s0++){
                for(unsigned int s1 = 0;s1 < 3;s1++){
                    for(unsigned int s2 = 0;s2 < 3;s2++){
                        unsigned int rightSubTree = 0;
                        unsigned int right = 0;
                        unsigned int perm = 0;
                        if(signs[s2]){//for (2) subtree
                            right = 0b1110;// 3 -> 2
                        } else{
                            right = 0b1011;// 2 -> 3
                        }

                        if(signs[s1]){//for (1 + 2) subtree (1 is left)
                            rightSubTree = (right << 2) | (0b01);// (2,3) -> 1
                        } else{
                            rightSubTree = ((0b01) << 4) | right;// 1 -> (2,3)
                        }

                        if(signs[s0]){
                            perm = (rightSubTree << 2) | (0b00);
                        } else{
                            perm = (0 << 6) + rightSubTree;
                        }
                        tempLUT[rs][s0 + s1 * 3 + s2 * 9 + 1 * 27] = PermToIndex(perm);
                    }
                }
            }


            //unsymetric topology 2  -> (1 + (2 + 1))
            for(unsigned int s0 = 0;s0 < 3;s0++){
                for(unsigned int s1 = 0;s1 < 3;s1++){
                    for(unsigned int s2 = 0;s2 < 3;s2++){
                        unsigned int SubTree = 0;
                        unsigned int left = 0;
                        unsigned int perm = 0;
                        if(signs[s2]){//for (2) subtree
                            left = 0b1001;// 2 -> 1
                        } else{
                            left = 0b0110;// 1 -> 2
                        }

                        if(signs[s1]){//for (1 + 2) subtree (1 is left)
                            SubTree = ((0b11) << 4) | left;// (2,3) -> 1
                        } else{
                            SubTree = (left << 2) | (0b11);// 1 -> (2,3)
                        }

                        if(signs[s0]){
                            perm = (SubTree << 2) | (0b00);
                        } else{
                            perm = (0 << 6) + SubTree;
                        }


                        tempLUT[rs][s0 + s1 * 3 + s2 * 9 + 2 * 27] = PermToIndex(perm);
                    }
                }
            }

            //unsymetric topology 3 -> ((2 + 1) + 1)
            for(unsigned int s0 = 0;s0 < 3;s0++){
                for(unsigned int s1 = 0;s1 < 3;s1++){
                    for(unsigned int s2 = 0;s2 < 3;s2++){
                        unsigned int SubTree = 0;
                        unsigned int left = 0;
                        unsigned int perm = 0;
                        if(signs[s2]){
                            left = 0b0100;
                        } else{
                            left = 0b0001;
                        }

                        if(signs[s1]){
                            SubTree = ((0b10) << 4) | (left);
                        } else{
                            SubTree = (left << 2) | (0b10);
                        }

                        if(signs[s0]){
                            perm = ((0b11) << 6) | (SubTree);
                        } else{
                            perm = (SubTree << 2) | (0b11);
                        }

                        tempLUT[rs][s0 + s1 * 3 + s2 * 9 + 3 * 27] = PermToIndex(perm);
                    }
                }
            }

            //unsymetric topology 4 -> ((1 + 2) + 1)
            for(unsigned int s0 = 0;s0 < 3;s0++){
                for(unsigned int s1 = 0;s1 < 3;s1++){
                    for(unsigned int s2 = 0;s2 < 3;s2++){
                        unsigned int SubTree = 0;
                        unsigned int left = 0;
                        unsigned int perm = 0;
                        if(signs[s2]){
                            left = 0b1001;
                        } else{
                            left = 0b0110;
                        }

                        if(signs[s1]){
                            SubTree = (left << 2) | (0b00);
                        } else{
                            SubTree = ((0b00) << 4) | left;
                        }

                        if(signs[s0]){
                            perm = ((0b11) << 6) | (SubTree);
                        } else{
                            perm = (SubTree << 2) | (0b11);
                        }

                        tempLUT[rs][s0 + s1 * 3 + s2 * 9 + 4 * 27] = PermToIndex(perm);
                    }
                }
            }
        }
        return tempLUT;
        }();
    static constexpr std::array<std::array<unsigned char, 24>, 16> maskLUT = [](){
        std::array<std::array<unsigned char, 24>, 16> tempMaskLUT;
        for(unsigned int mask = 0;mask < 16;mask++){
            for(unsigned int perm = 0;perm < 24;perm++){
                const unsigned int order = PermToIndexLUT[perm];
                unsigned char ans = 0;
                for(int n = 6;n >= 0;n -= 2){
                    const unsigned char idx = (order >> n) & 0x3;
                    if(mask & (1 << idx)){
                        //0321
                        //stack = 1 2 3 -> 3 is first to intersect
                        ans <<= 2;
                        ans |= idx;
                    }
                }
                tempMaskLUT[mask][perm] = ans;
            }
        }
        return tempMaskLUT;
        }();
public:
    virtual ~BVH4() = default;
    BVH4() = default;

    BVH4(const std::vector<T>& prims){
        auto start = std::chrono::high_resolution_clock::now();
        std::cout << "Building BVH4" << std::endl;

        nodes.reserve(prims.size() * 2 - 1);
        BVHBase<T>::primitives.reserve(prims.size());

        std::vector<PrimitiveInfo> primitiveInfo;
        primitiveInfo.reserve(prims.size());

        for(std::size_t i = 0;i < prims.size();i++){
            if constexpr(requires { prims[i]->BoundingBox(); }){
                AABB bbox = prims[i]->BoundingBox();
                primitiveInfo.emplace_back(i, bbox);
            } else{
                AABB bbox = prims[i].BoundingBox();
                primitiveInfo.emplace_back(i, bbox);
            }
        }

        //buildBVHBase accepty num of tris in leaf?
        //if std::is_base_of<TrianglePrimitive> 
        //then we gather 8 shapes with GetShape
        //we then put all the data in the float arrays
        //pass them to Intersect8 and get result from it
        std::vector<BVH_NODE> BVHnodes = BVHBase<T>::BuildBVHBase(primitiveInfo);
        if(BVHnodes.empty()){
            BVHBase<T>::BVHbbox = AABB {};
        } else{
            BVHBase<T>::BVHbbox = BVHnodes[0].bbox;
        }


        rootNode = buildBVH4(BVHnodes, 0);
        std::cout << BVHnodes.size() << "  " << nodes.size() << "....\n";

        this->nodes.shrink_to_fit();
        for(const auto& info : primitiveInfo){
            BVHBase<T>::primitives.emplace_back(prims[info.index]);
        }

        auto duration = std::chrono::high_resolution_clock::now() - start;
        std::cout << "BVH4 built in: " << std::chrono::duration_cast<std::chrono::milliseconds>(duration) << std::endl;
    }

    BVH4_NODE buildBVH4(const std::vector<BVH_NODE>& BVHnodes, uint32_t nodeIdx = 0){
        const BVH_NODE& n = BVHnodes[nodeIdx];
        unsigned char active = 0;
        unsigned char perm = 0;
        if(n.isLeaf()){
            return BVH4_NODE { static_cast<unsigned char>(n.count),0b0000,0,n.right };
        }

        uint32_t index = nodes.size();
        BVH4_CLUSTER& node = nodes.emplace_back();
        const BVH_NODE& left = BVHnodes[n.left];
        const BVH_NODE& right = BVHnodes[n.right];

        if(left.isLeaf() && right.isLeaf()){
            /*
                     n
                   /   \
               node0   node1
            */
            BVH4_NODE node0 = buildBVH4(BVHnodes, n.left);
            BVH4_NODE node2 = buildBVH4(BVHnodes, n.right);

            node.SetDataAtIndex(left.bbox, 0);
            node.SetDataAtIndex(right.bbox, 2);
            node.children[0] = node0;
            node.children[2] = node2;
            active = 0b0101;
            perm = n.axis + 0 * 3 + 0 * 9 + 0 * 27;
        } else if(left.isLeaf()){
            BVH_NODE lowerLeft = BVHnodes[right.left];
            BVH_NODE lowerRight = BVHnodes[right.right];
            if(lowerLeft.isLeaf() && lowerRight.isLeaf()){
                /*
                         n
                       /  \
                   node0    right
                           /    \
                       node1    node2
                */
                BVH4_NODE node0 = buildBVH4(BVHnodes, n.left);
                BVH4_NODE node2 = buildBVH4(BVHnodes, right.left);
                BVH4_NODE node3 = buildBVH4(BVHnodes, right.right);

                node.SetDataAtIndex(left.bbox, 0);
                node.SetDataAtIndex(lowerLeft.bbox, 2);
                node.SetDataAtIndex(lowerRight.bbox, 3);

                node.children[0] = node0;
                node.children[2] = node2;
                node.children[3] = node3;

                active = 0b1101;
                perm = n.axis + 0 * 3 + right.axis * 9 + 0 * 27;
            } else if(lowerLeft.isLeaf()){
                /*
                         n
                       /  \
                   node0    right
                           /    \
                       node1    lowerRight
                                /        \
                            node2        node3
                */
                BVH_NODE l = BVHnodes[lowerRight.left];
                BVH_NODE r = BVHnodes[lowerRight.right];
                BVH4_NODE node0 = buildBVH4(BVHnodes, n.left);
                BVH4_NODE node1 = buildBVH4(BVHnodes, right.left);
                BVH4_NODE node2 = buildBVH4(BVHnodes, lowerRight.left);
                BVH4_NODE node3 = buildBVH4(BVHnodes, lowerRight.right);

                node.SetDataAtIndex(left.bbox, 0);
                node.SetDataAtIndex(lowerLeft.bbox, 1);
                node.SetDataAtIndex(l.bbox, 2);
                node.SetDataAtIndex(r.bbox, 3);

                node.children[0] = node0;
                node.children[1] = node1;
                node.children[2] = node2;
                node.children[3] = node3;

                active = 0b1111;
                perm = n.axis + right.axis * 3 + lowerRight.axis * 9 + 1 * 27;
            } else{
                /*
                         n
                       /  \
                   node0    right
                           /     \
                     lowerLeft   node4
                      /     \
                  node1     node2
                */
                BVH_NODE l = BVHnodes[lowerLeft.left];
                BVH_NODE r = BVHnodes[lowerLeft.right];
                BVH4_NODE node0 = buildBVH4(BVHnodes, n.left);
                BVH4_NODE node1 = buildBVH4(BVHnodes, lowerLeft.left);
                BVH4_NODE node2 = buildBVH4(BVHnodes, lowerLeft.right);
                BVH4_NODE node3 = buildBVH4(BVHnodes, right.right);


                node.SetDataAtIndex(left.bbox, 0);
                node.SetDataAtIndex(l.bbox, 1);
                node.SetDataAtIndex(r.bbox, 2);
                node.SetDataAtIndex(lowerRight.bbox, 3);

                node.children[0] = node0;
                node.children[1] = node1;
                node.children[2] = node2;
                node.children[3] = node3;

                active = 0b1111;
                perm = n.axis + right.axis * 3 + lowerLeft.axis * 9 + 2 * 27;
            }
        } else if(right.isLeaf()){
            BVH_NODE lowerLeft = BVHnodes[left.left];
            BVH_NODE lowerRight = BVHnodes[left.right];
            if(lowerLeft.isLeaf() && lowerRight.isLeaf()){
                /*
                              n
                           /    \
                       left      node2
                       /   \
                  node0  node1
                */
                BVH4_NODE node0 = buildBVH4(BVHnodes, left.left);
                BVH4_NODE node1 = buildBVH4(BVHnodes, left.right);
                BVH4_NODE node2 = buildBVH4(BVHnodes, n.right);

                node.SetDataAtIndex(lowerLeft.bbox, 0);
                node.SetDataAtIndex(lowerRight.bbox, 1);
                node.SetDataAtIndex(right.bbox, 2);

                node.children[0] = node0;
                node.children[1] = node1;
                node.children[2] = node2;

                active = 0b0111;
                perm = n.axis + left.axis * 3 + 0 * 9 + 0 * 27;
            } else if(lowerLeft.isLeaf()){
                /*
                                n
                             /    \
                         left      node3
                        /    \
                   node0  lowerRight
                           /     \
                       node1     node2
                */
                BVH_NODE l = BVHnodes[lowerRight.left];
                BVH_NODE r = BVHnodes[lowerRight.right];
                BVH4_NODE node0 = buildBVH4(BVHnodes, left.left);
                BVH4_NODE node1 = buildBVH4(BVHnodes, lowerRight.left);
                BVH4_NODE node2 = buildBVH4(BVHnodes, lowerRight.right);
                BVH4_NODE node3 = buildBVH4(BVHnodes, n.right);

                node.SetDataAtIndex(lowerLeft.bbox, 0);
                node.SetDataAtIndex(l.bbox, 1);
                node.SetDataAtIndex(r.bbox, 2);
                node.SetDataAtIndex(right.bbox, 3);

                node.children[0] = node0;
                node.children[1] = node1;
                node.children[2] = node2;
                node.children[3] = node3;

                active = 0b1111;
                perm = n.axis + left.axis * 3 + lowerRight.axis * 9 + 4 * 27;
            } else{
                /*
                                   n
                                /    \
                            left      node3
                           /    \
                      lowerLeft  node2
                       /     \
                  node0     node1
                */
                BVH_NODE l = BVHnodes[lowerLeft.left];
                BVH_NODE r = BVHnodes[lowerLeft.right];
                BVH4_NODE node0 = buildBVH4(BVHnodes, lowerLeft.left);
                BVH4_NODE node1 = buildBVH4(BVHnodes, lowerLeft.right);
                BVH4_NODE node2 = buildBVH4(BVHnodes, left.right);
                BVH4_NODE node3 = buildBVH4(BVHnodes, n.right);

                node.SetDataAtIndex(l.bbox, 0);
                node.SetDataAtIndex(r.bbox, 1);
                node.SetDataAtIndex(lowerRight.bbox, 2);
                node.SetDataAtIndex(right.bbox, 3);

                node.children[0] = node0;
                node.children[1] = node1;
                node.children[2] = node2;
                node.children[3] = node3;

                active = 0b1111;
                perm = n.axis + left.axis * 3 + lowerLeft.axis * 9 + 3 * 27;
            }
        } else{
            /*
                           n
                        /     \
                    left      right
                   /    \     /    \
               node0  node1 node2 node3
            */
            BVH_NODE A = BVHnodes[left.left];
            BVH_NODE B = BVHnodes[left.right];
            BVH_NODE C = BVHnodes[right.left];
            BVH_NODE D = BVHnodes[right.right];

            BVH4_NODE node0 = buildBVH4(BVHnodes, left.left);
            BVH4_NODE node1 = buildBVH4(BVHnodes, left.right);
            BVH4_NODE node2 = buildBVH4(BVHnodes, right.left);
            BVH4_NODE node3 = buildBVH4(BVHnodes, right.right);

            node.SetDataAtIndex(A.bbox, 0);
            node.SetDataAtIndex(B.bbox, 1);
            node.SetDataAtIndex(C.bbox, 2);
            node.SetDataAtIndex(D.bbox, 3);

            node.children[0] = node0;
            node.children[1] = node1;
            node.children[2] = node2;
            node.children[3] = node3;

            active = 0b1111;
            perm = n.axis + left.axis * 3 + right.axis * 9 + 0 * 27;
        }
        return BVH4_NODE { static_cast<unsigned char>(n.count),active,perm,index };
    }

    bool IntersectPred(const Ray& ray, float max = std::numeric_limits<float>::infinity()) const final{
        static const __m128 Eps = _mm_set1_ps(shadowEpsilon);
        //origin
        const __m128 Ox4 = _mm_set1_ps(ray.origin.x);
        const __m128 Oy4 = _mm_set1_ps(ray.origin.y);
        const __m128 Oz4 = _mm_set1_ps(ray.origin.z);

        //inverse of dir
        const __m128 rDx4 = _mm_set1_ps(ray.inv_dir.x);
        const __m128 rDy4 = _mm_set1_ps(ray.inv_dir.y);
        const __m128 rDz4 = _mm_set1_ps(ray.inv_dir.z);

        //ray signs
        //lowest bit is x, becouse we use axis0 = x, see building of LUT
        //const unsigned char signs = ((ray.dir[2] < 0) << 2) | ((ray.dir[1] < 0) << 1) | (ray.dir[0] < 0);
        const __m128 maxt = _mm_set1_ps(max);
        BVH4_NODE stack[32];
        int sp = 0;
        stack[sp++] = rootNode;
        while(sp){
            BVH4_NODE bvh4node = stack[--sp];

            if(bvh4node.active != 0){
                const BVH4_CLUSTER& node = nodes[bvh4node.ClusterIdx];
                //not leaf

                //load everything into simd
                const __m128 xmin4 = _mm_load_ps(node.xmin);
                const __m128 xmax4 = _mm_load_ps(node.xmax);
                const __m128 ymin4 = _mm_load_ps(node.ymin);
                const __m128 ymax4 = _mm_load_ps(node.ymax);
                const __m128 zmin4 = _mm_load_ps(node.zmin);
                const __m128 zmax4 = _mm_load_ps(node.zmax);

                //const __m128 t1 = _mm_mul_ps(_mm_sub_ps(bmin4,ray.O4),ray.rD4 );
                const __m128 tx1 = _mm_mul_ps(_mm_sub_ps(xmin4, Ox4), rDx4);
                const __m128 ty1 = _mm_mul_ps(_mm_sub_ps(ymin4, Oy4), rDy4);
                const __m128 tz1 = _mm_mul_ps(_mm_sub_ps(zmin4, Oz4), rDz4);

                //const __m128 t2 = _mm_mul_ps(_mm_sub_ps(bmax4,ray.O4),ray.rD4 );
                const __m128 tx2 = _mm_mul_ps(_mm_sub_ps(xmax4, Ox4), rDx4);
                const __m128 ty2 = _mm_mul_ps(_mm_sub_ps(ymax4, Oy4), rDy4);
                const __m128 tz2 = _mm_mul_ps(_mm_sub_ps(zmax4, Oz4), rDz4);

                //const __m128 vmax4 = _mm_max_ps(t1,t2);
                const __m128 tmaxx4 = _mm_max_ps(tx1, tx2);
                const __m128 tmaxy4 = _mm_max_ps(ty1, ty2);
                const __m128 tmaxz4 = _mm_max_ps(tz1, tz2);

                //const __m128 vmin4 = _mm_min_ps(t1,t2);
                const __m128 tminx4 = _mm_min_ps(tx1, tx2);
                const __m128 tminy4 = _mm_min_ps(ty1, ty2);
                const __m128 tminz4 = _mm_min_ps(tz1, tz2);

                const __m128 tminxy = _mm_max_ps(tminx4, tminy4);
                const __m128 tEntry = _mm_max_ps(tminxy, tminz4);//tEntry

                const __m128 tmaxxy = _mm_min_ps(tmaxx4, tmaxy4);
                const __m128 tExit = _mm_min_ps(tmaxxy, tmaxz4);//tExit



                /*
                __m128i active_i = _mm_setr_epi32(
                    (node.active & 1) ? -1 : 0,               // lane 0 (x)
                    ((node.active >> 1) & 1) ? -1 : 0,        // lane 1 (y)
                    ((node.active >> 2) & 1) ? -1 : 0,        // lane 2 (z)
                    ((node.active >> 3) & 1) ? -1 : 0         // lane 3
                );
                __m128 activeMask = _mm_castsi128_ps(active_i);
                */

                __m128 hitmask = _mm_cmpge_ps(tExit, Eps);//tExit >= 0
                hitmask = _mm_and_ps(hitmask, _mm_cmplt_ps(tEntry, maxt));//tEntry <= maxT
                hitmask = _mm_and_ps(hitmask, _mm_cmple_ps(tEntry, tExit));//tEntry <= tExit
                //hitmask = _mm_and_ps(hitmask,activeMask);

                unsigned int mask = _mm_movemask_ps(hitmask);


                for(int i = 0;i < 4;i++){
                    if(mask & 1)stack[sp++] = node.children[i];
                    mask >>= 1;
                }
                /*
                int bitcnt = std::popcount(mask);
                const unsigned char orderIdx = LUT[signs][bvh4node.perm];
                unsigned char order = maskLUT[mask][orderIdx];

                while(bitcnt--!=0){
                    stack[sp++]=node.children[order&0x3];
                    order>>=2;
                }*/

            } else if(BVHBase<T>::intersectPrimitivesPred(ray, max, bvh4node.ClusterIdx, bvh4node.count)){
                return true;
            }
        }
        return false;
    }

    bool Intersect(const Ray& ray, SurfaceInteraction& interaction, float max = std::numeric_limits<float>::infinity()) const final{
        static const __m128 Eps = _mm_set1_ps(shadowEpsilon);
        //origin
        const __m128 Ox4 = _mm_set1_ps(ray.origin.x);
        const __m128 Oy4 = _mm_set1_ps(ray.origin.y);
        const __m128 Oz4 = _mm_set1_ps(ray.origin.z);

        //inverse of dir
        const __m128 rDx4 = _mm_set1_ps(ray.inv_dir.x);
        const __m128 rDy4 = _mm_set1_ps(ray.inv_dir.y);
        const __m128 rDz4 = _mm_set1_ps(ray.inv_dir.z);

        //ray signs
        //lowest bit is x, becouse we use axis0 = x, see building of LUT
        const unsigned char signs = ((ray.dir[2] < 0) << 2) | ((ray.dir[1] < 0) << 1) | (ray.dir[0] < 0);
        __m128 maxt = _mm_set1_ps(max);
        BVH4_NODE stack[32];
        float entryDist[32];
        int sp = 0;
        entryDist[sp] = 0;
        stack[sp++] = rootNode;
        bool hit = false;
        while(sp){
            if(entryDist[--sp] > max)continue;
            BVH4_NODE bvh4node = stack[sp];
            if(bvh4node.active != 0){
                const BVH4_CLUSTER& node = nodes[bvh4node.ClusterIdx];
                //not leaf

                //load everything into simd
                const __m128 xmin4 = _mm_load_ps(node.xmin);
                const __m128 xmax4 = _mm_load_ps(node.xmax);
                const __m128 ymin4 = _mm_load_ps(node.ymin);
                const __m128 ymax4 = _mm_load_ps(node.ymax);
                const __m128 zmin4 = _mm_load_ps(node.zmin);
                const __m128 zmax4 = _mm_load_ps(node.zmax);

                //const __m128 t1 = _mm_mul_ps(_mm_sub_ps(bmin4,ray.O4),ray.rD4 );
                const __m128 tx1 = _mm_mul_ps(_mm_sub_ps(xmin4, Ox4), rDx4);
                const __m128 ty1 = _mm_mul_ps(_mm_sub_ps(ymin4, Oy4), rDy4);
                const __m128 tz1 = _mm_mul_ps(_mm_sub_ps(zmin4, Oz4), rDz4);

                //const __m128 t2 = _mm_mul_ps(_mm_sub_ps(bmax4,ray.O4),ray.rD4 );
                const __m128 tx2 = _mm_mul_ps(_mm_sub_ps(xmax4, Ox4), rDx4);
                const __m128 ty2 = _mm_mul_ps(_mm_sub_ps(ymax4, Oy4), rDy4);
                const __m128 tz2 = _mm_mul_ps(_mm_sub_ps(zmax4, Oz4), rDz4);

                //const __m128 vmax4 = _mm_max_ps(t1,t2);
                const __m128 tmaxx4 = _mm_max_ps(tx1, tx2);
                const __m128 tmaxy4 = _mm_max_ps(ty1, ty2);
                const __m128 tmaxz4 = _mm_max_ps(tz1, tz2);

                //const __m128 vmin4 = _mm_min_ps(t1,t2);
                const __m128 tminx4 = _mm_min_ps(tx1, tx2);
                const __m128 tminy4 = _mm_min_ps(ty1, ty2);
                const __m128 tminz4 = _mm_min_ps(tz1, tz2);

                const __m128 tminxy = _mm_max_ps(tminx4, tminy4);
                const __m128 tEntry = _mm_max_ps(tminxy, tminz4);//tEntry

                const __m128 tmaxxy = _mm_min_ps(tmaxx4, tmaxy4);
                const __m128 tExit = _mm_min_ps(tmaxxy, tmaxz4);//tExit



                /*
                __m128i active_i = _mm_setr_epi32(
                    (node.active & 1) ? -1 : 0,               // lane 0 (x)
                    ((node.active >> 1) & 1) ? -1 : 0,        // lane 1 (y)
                    ((node.active >> 2) & 1) ? -1 : 0,        // lane 2 (z)
                    ((node.active >> 3) & 1) ? -1 : 0         // lane 3
                );
                __m128 activeMask = _mm_castsi128_ps(active_i);
                */

                __m128 hitmask = _mm_cmpge_ps(tExit, Eps);//tExit >= 0
                hitmask = _mm_and_ps(hitmask, _mm_cmplt_ps(tEntry, maxt));//tEntry <= maxT
                hitmask = _mm_and_ps(hitmask, _mm_cmple_ps(tEntry, tExit));//tEntry <= tExit
                //hitmask = _mm_and_ps(hitmask,activeMask);

                alignas(16) float distances[4];
                _mm_store_ps(distances, tEntry);

                const unsigned int mask = _mm_movemask_ps(hitmask);
                int bitcnt = std::popcount(mask);
                const unsigned char orderIdx = LUT[signs][bvh4node.perm];
                unsigned char order = maskLUT[mask][orderIdx];

                while(bitcnt-- != 0){
                    entryDist[sp] = distances[order & 0x3];
                    stack[sp++] = node.children[order & 0x3];
                    order >>= 2;
                }
            } else{
                hit |= BVHBase<T>::intersectPrimitives(ray, max, bvh4node.ClusterIdx, bvh4node.count, interaction);
                maxt = _mm_set1_ps(max);
            }
        }
        return hit;
    }


protected:
    BVH4_NODE rootNode;
    std::vector<BVH4_CLUSTER> nodes;
};

using BLAS4 = BVH4<GeometricPrimitive>;
using TLAS4 = BVH4<std::shared_ptr<Primitive>>;
#endif


//during BVH8 construction keep all nodes (even the inner temporary ones)
//after then find order for each ray sign combination
//place the final nodes into tbe BVH8_CLUSTER

struct BVH8_NODE{
    uint32_t perm = 0;//top 8 bytes hold count -> if count != 0 -> leaf
    //1. way
    //(perm&0x7) << (3* ray sign) to find position of that child
    //2. way 
    // child[ray sign].perm holds all rodering
    // this could have been used in BVH4 becouse we need 8 perm vectors -> aka >=8 children
    uint32_t ClusterIdx = 0;//if leaf then this holds triangle index
};
//if we use 1. we can just >> (3*ray sign(0bXXX)) to get the nodes index in bottom 3 bits
//we can maybe do this with vector instruction (rotate)
// to do this we need AVX 512

struct alignas(256) BVH8_CLUSTER{
    float xmin[8] { 0 }, xmax[8] { 0 };
    float ymin[8] { 0 }, ymax[8] { 0 };
    float zmin[8] { 0 }, zmax[8] { 0 };

    BVH8_NODE children[8] {};

    void SetDataAtIndex(const AABB& bbox, int index){
        xmin[index] = bbox.min[0];
        ymin[index] = bbox.min[1];
        zmin[index] = bbox.min[2];
        xmax[index] = bbox.max[0];
        ymax[index] = bbox.max[1];
        zmax[index] = bbox.max[2];
    }
};

#if defined(__AVX__) || defined(_M_AMD64) || defined(_M_X64)
template <typename T>
class BVH8 : public BVHBase<T>{
    using typename BVHBase<T>::PrimitiveInfo;
public:
    virtual ~BVH8() = default;
    BVH8() = default;

    BVH8(const std::vector<T>& prims){
        auto start = std::chrono::high_resolution_clock::now();
        std::cout << "Building BVH8" << std::endl;

        nodes.reserve(prims.size() * 2 - 1);
        BVHBase<T>::primitives.reserve(prims.size());

        std::vector<PrimitiveInfo> primitiveInfo;
        primitiveInfo.reserve(prims.size());

        for(std::size_t i = 0;i < prims.size();i++){
            if constexpr(requires { prims[i]->BoundingBox(); }){
                AABB bbox = prims[i]->BoundingBox();
                primitiveInfo.emplace_back(i, bbox);
            } else{
                AABB bbox = prims[i].BoundingBox();
                primitiveInfo.emplace_back(i, bbox);
            }
        }

        //buildBVHBase accepty num of tris in leaf?
        //if std::is_base_of<TrianglePrimitive> 
        //then we gather 8 shapes with GetShape
        //we then put all the data in the float arrays
        //pass them to Intersect8 and get result from it
        std::vector<BVH_NODE> BVHnodes = BVHBase<T>::BuildBVHBase(primitiveInfo);
        if(BVHnodes.empty()){
            BVHBase<T>::BVHbbox = AABB {};
        } else{
            BVHBase<T>::BVHbbox = BVHnodes[0].bbox;
        }


        rootNode = buildBVH8(BVHnodes, 0);
        std::cout << BVHnodes.size() << "  " << nodes.size() << "....\n";

        this->nodes.shrink_to_fit();
        for(const auto& info : primitiveInfo){
            BVHBase<T>::primitives.emplace_back(prims[info.index]);
        }

        auto duration = std::chrono::high_resolution_clock::now() - start;
        std::cout << "BVH8 built in: " << std::chrono::duration_cast<std::chrono::milliseconds>(duration) << std::endl;
    }

    BVH8_NODE buildBVH8(const std::vector<BVH_NODE>& BVHnodes, uint32_t nodeIdx = 0){
        const BVH_NODE& n = BVHnodes[nodeIdx];
        if(n.isLeaf()){
            return BVH8_NODE { static_cast<uint32_t>(n.count) << 24,n.right };
        }

        uint32_t index = nodes.size();
        BVH8_CLUSTER& node = nodes.emplace_back();

        std::array<BVH_NODE, 8> leaf_nodes = { n };
        std::array<uint32_t, 8> leaf_indices = { nodeIdx };
        int end = 1;
        while(end < 8){
            int k = -1;
            float maxArea = -1;
            for(int i = 0;i < end;i++){
                if(!leaf_nodes[i].isLeaf() && leaf_nodes[i].bbox.Area() >= maxArea){
                    k = i;
                    maxArea = leaf_nodes[i].bbox.Area();
                }
            }
            if(k == -1)break;
            leaf_nodes[end] = BVHnodes[leaf_nodes[k].right];
            leaf_indices[end] = leaf_nodes[k].right;
            leaf_indices[k] = leaf_nodes[k].left;
            leaf_nodes[k] = BVHnodes[leaf_nodes[k].left];
            ++end;
        }


        for(int i = 0;i < end;i++){
            node.SetDataAtIndex(leaf_nodes[i].bbox, i);
            node.children[i] = buildBVH8(BVHnodes, leaf_indices[i]);
        }

        for(unsigned int rs = 0;rs < 8;rs++){
            const unsigned int signs[3] = { rs & 1 , (rs >> 1) & 1, (rs >> 2) };
            uint32_t order = 0xFFFFFF;
            buildPerm(signs, leaf_indices, end, nodeIdx, BVHnodes, order);
            node.children[rs].perm |= order & 0xFFFFFF;
        }

        return BVH8_NODE { 0,index };
    }

    void buildPerm(const unsigned int signs[3], const std::array<uint32_t, 8>& leaf_indices, int end, uint32_t curr, const std::vector<BVH_NODE>& BVHnodes, uint32_t& order) const{
        for(uint32_t i = 0;i < end;i++){
            if(leaf_indices[i] == curr){
                order <<= 3;
                order |= i;
                return;
            }
        }
        const BVH_NODE& n = BVHnodes[curr];
        if(signs[n.axis]){
            buildPerm(signs, leaf_indices, end, n.right, BVHnodes, order);
            buildPerm(signs, leaf_indices, end, n.left, BVHnodes, order);
        } else{
            buildPerm(signs, leaf_indices, end, n.left, BVHnodes, order);
            buildPerm(signs, leaf_indices, end, n.right, BVHnodes, order);
        }
    }

    bool IntersectPred(const Ray& ray, float max = std::numeric_limits<float>::infinity()) const final{
        static const __m256 Eps = _mm256_set1_ps(shadowEpsilon);
        //origin
        const __m256 Ox4 = _mm256_set1_ps(ray.origin.x);
        const __m256 Oy4 = _mm256_set1_ps(ray.origin.y);
        const __m256 Oz4 = _mm256_set1_ps(ray.origin.z);

        //inverse of dir
        const __m256 rDx4 = _mm256_set1_ps(ray.inv_dir.x);
        const __m256 rDy4 = _mm256_set1_ps(ray.inv_dir.y);
        const __m256 rDz4 = _mm256_set1_ps(ray.inv_dir.z);

        //ray signs
        //lowest bit is x, becouse we use axis0 = x, see building of LUT
        //const unsigned char signs = ((ray.dir[2] < 0) << 2) | ((ray.dir[1] < 0) << 1) | (ray.dir[0] < 0);
        const __m256 maxt = _mm256_set1_ps(max);
        BVH8_NODE stack[64];
        int sp = 0;
        stack[sp++] = rootNode;
        while(sp){
            BVH8_NODE bvh8node = stack[--sp];

            if((bvh8node.perm >> 24) == 0){
                const BVH8_CLUSTER& node = nodes[bvh8node.ClusterIdx];
                //not leaf

                //load everything into simd
                const __m256 xmin4 = _mm256_load_ps(node.xmin);
                const __m256 xmax4 = _mm256_load_ps(node.xmax);
                const __m256 ymin4 = _mm256_load_ps(node.ymin);
                const __m256 ymax4 = _mm256_load_ps(node.ymax);
                const __m256 zmin4 = _mm256_load_ps(node.zmin);
                const __m256 zmax4 = _mm256_load_ps(node.zmax);

                //const __m256 t1 = _mm256_mul_ps(_mm256_sub_ps(bmin4,ray.O4),ray.rD4 );
                const __m256 tx1 = _mm256_mul_ps(_mm256_sub_ps(xmin4, Ox4), rDx4);
                const __m256 ty1 = _mm256_mul_ps(_mm256_sub_ps(ymin4, Oy4), rDy4);
                const __m256 tz1 = _mm256_mul_ps(_mm256_sub_ps(zmin4, Oz4), rDz4);

                //const __m256 t2 = _mm256_mul_ps(_mm256_sub_ps(bmax4,ray.O4),ray.rD4 );
                const __m256 tx2 = _mm256_mul_ps(_mm256_sub_ps(xmax4, Ox4), rDx4);
                const __m256 ty2 = _mm256_mul_ps(_mm256_sub_ps(ymax4, Oy4), rDy4);
                const __m256 tz2 = _mm256_mul_ps(_mm256_sub_ps(zmax4, Oz4), rDz4);

                //const __m256 vmax4 = _mm256_max_ps(t1,t2);
                const __m256 tmaxx4 = _mm256_max_ps(tx1, tx2);
                const __m256 tmaxy4 = _mm256_max_ps(ty1, ty2);
                const __m256 tmaxz4 = _mm256_max_ps(tz1, tz2);

                //const __m256 vmin4 = _mm256_min_ps(t1,t2);
                const __m256 tminx4 = _mm256_min_ps(tx1, tx2);
                const __m256 tminy4 = _mm256_min_ps(ty1, ty2);
                const __m256 tminz4 = _mm256_min_ps(tz1, tz2);

                const __m256 tminxy = _mm256_max_ps(tminx4, tminy4);
                const __m256 tEntry = _mm256_max_ps(tminxy, tminz4);//tEntry

                const __m256 tmaxxy = _mm256_min_ps(tmaxx4, tmaxy4);
                const __m256 tExit = _mm256_min_ps(tmaxxy, tmaxz4);//tExit



                __m256 hitmask = _mm256_cmp_ps(tExit, Eps, _CMP_GE_OQ);//tExit >= 0
                hitmask = _mm256_and_ps(hitmask, _mm256_cmp_ps(tEntry, maxt, _CMP_LT_OQ));//tEntry <= maxT
                hitmask = _mm256_and_ps(hitmask, _mm256_cmp_ps(tEntry, tExit, _CMP_LE_OQ));//tEntry <= tExit

                unsigned int mask = _mm256_movemask_ps(hitmask);


                for(int i = 0;i < 8;i++){
                    if((mask >> i) & 1)stack[sp++] = node.children[i];
                }

            } else if(BVHBase<T>::intersectPrimitivesPred(ray, max, bvh8node.ClusterIdx, bvh8node.perm >> 24)){
                return true;
            }
        }
        return false;
    }

    bool Intersect(const Ray& ray, SurfaceInteraction& interaction, float max = std::numeric_limits<float>::infinity()) const final{
        static const __m256 Eps = _mm256_set1_ps(shadowEpsilon);
        //origin
        const __m256 Ox8 = _mm256_set1_ps(ray.origin.x);
        const __m256 Oy8 = _mm256_set1_ps(ray.origin.y);
        const __m256 Oz8 = _mm256_set1_ps(ray.origin.z);

        //inverse of dir
        const __m256 rDx8 = _mm256_set1_ps(ray.inv_dir.x);
        const __m256 rDy8 = _mm256_set1_ps(ray.inv_dir.y);
        const __m256 rDz8 = _mm256_set1_ps(ray.inv_dir.z);

        //ray signs
        //lowest bit is x, becouse we use axis0 = x, see building of LUT
        const unsigned char signs = ((ray.dir[2] < 0) << 2) | ((ray.dir[1] < 0) << 1) | (ray.dir[0] < 0);
        __m256 maxt = _mm256_set1_ps(max);
        BVH8_NODE stack[64];
        float entryDist[64];
        int sp = 0;
        entryDist[sp] = 0;
        stack[sp++] = rootNode;
        bool hit = false;
        while(sp){
            if(entryDist[--sp] > max)continue;
            BVH8_NODE bvh8node = stack[sp];
            if((bvh8node.perm >> 24) == 0){
                const BVH8_CLUSTER& node = nodes[bvh8node.ClusterIdx];
                //not leaf

                //load everything into simd
                const __m256 xmin8 = _mm256_load_ps(node.xmin);
                const __m256 xmax8 = _mm256_load_ps(node.xmax);
                const __m256 ymin8 = _mm256_load_ps(node.ymin);
                const __m256 ymax8 = _mm256_load_ps(node.ymax);
                const __m256 zmin8 = _mm256_load_ps(node.zmin);
                const __m256 zmax8 = _mm256_load_ps(node.zmax);
                uint32_t order = node.children[signs].perm;

                //const __m256 t1 = _mm256_mul_ps(_mm256_sub_ps(bmin8,ray.O8),ray.rD8 );
                const __m256 tx1 = _mm256_mul_ps(_mm256_sub_ps(xmin8, Ox8), rDx8);
                const __m256 ty1 = _mm256_mul_ps(_mm256_sub_ps(ymin8, Oy8), rDy8);
                const __m256 tz1 = _mm256_mul_ps(_mm256_sub_ps(zmin8, Oz8), rDz8);

                //const __m256 t2 = _mm256_mul_ps(_mm256_sub_ps(bmax8,ray.O8),ray.rD8 );
                const __m256 tx2 = _mm256_mul_ps(_mm256_sub_ps(xmax8, Ox8), rDx8);
                const __m256 ty2 = _mm256_mul_ps(_mm256_sub_ps(ymax8, Oy8), rDy8);
                const __m256 tz2 = _mm256_mul_ps(_mm256_sub_ps(zmax8, Oz8), rDz8);

                //const __m256 vmax8 = _mm256_max_ps(t1,t2);
                const __m256 tmaxx8 = _mm256_max_ps(tx1, tx2);
                const __m256 tmaxy8 = _mm256_max_ps(ty1, ty2);
                const __m256 tmaxz8 = _mm256_max_ps(tz1, tz2);

                //const __m256 vmin8 = _mm256_min_ps(t1,t2);
                const __m256 tminx8 = _mm256_min_ps(tx1, tx2);
                const __m256 tminy8 = _mm256_min_ps(ty1, ty2);
                const __m256 tminz8 = _mm256_min_ps(tz1, tz2);

                const __m256 tminxy = _mm256_max_ps(tminx8, tminy8);
                const __m256 tEntry = _mm256_max_ps(tminxy, tminz8);//tEntry

                const __m256 tmaxxy = _mm256_min_ps(tmaxx8, tmaxy8);
                const __m256 tExit = _mm256_min_ps(tmaxxy, tmaxz8);//tExit


                __m256 hitmask = _mm256_cmp_ps(tExit, Eps, _CMP_GE_OQ);//tExit >= 0
                hitmask = _mm256_and_ps(hitmask, _mm256_cmp_ps(tEntry, maxt, _CMP_LT_OQ));//tEntry <= maxT -> maybe remove
                hitmask = _mm256_and_ps(hitmask, _mm256_cmp_ps(tEntry, tExit, _CMP_LE_OQ));//tEntry <= tExit
                //hitmask + permutation vector- -> _mm256_permute_ps ?
                alignas(32) float distances[8];
                _mm256_store_ps(distances, tEntry);

                unsigned int mask = _mm256_movemask_ps(hitmask);

                for(int i = 0;i < 8;i++){
                    unsigned char t = order & 0b111;
                    order >>= 3;
                    if((mask >> t) & 1){
                        entryDist[sp] = distances[t];
                        stack[sp++] = node.children[t];
                    }
                }
            } else{
                hit |= BVHBase<T>::intersectPrimitives(ray, max, bvh8node.ClusterIdx, bvh8node.perm >> 24, interaction);
                maxt = _mm256_set1_ps(max);
            }
        }
        return hit;
    }


protected:
    BVH8_NODE rootNode;
    std::vector<BVH8_CLUSTER> nodes;
};

using BLAS8 = BVH8<GeometricPrimitive>;
using TLAS8 = BVH8<std::shared_ptr<Primitive>>;
#endif