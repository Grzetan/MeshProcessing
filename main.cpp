#include <iostream>
#include <algorithm>
#include <math.h>
#include <map>
#include <array>
#include "list.h"
#include "happly.h"

typedef double real;
typedef std::array<real, 3> point3d;
typedef std::vector<size_t> polygonIndexes;
typedef std::vector<point3d> pointCloud3d;
typedef std::tuple<int, int, int> mapKey;
typedef std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, int>>> sharedMap;
typedef sharedMap::iterator mapIter;

point3d operator-(const point3d &p) {
    return point3d{-p[0], -p[1], -p[2]};
}

point3d CalculateCenterPointOfTriangle(const point3d &a, const point3d &b, const point3d &c) {
    point3d center = {(a[0] + b[0] + c[0]) / 3,
                      (a[1] + b[1] + c[1]) / 3,
                      (a[2] + b[2] + c[2]) / 3};
    return center;
}

real dotProduct(const point3d &a, const point3d &b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

point3d crossProduct(const point3d &a, const point3d &b) {
    point3d c = {0, 0, 0};
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = -(a[0] * b[2] - a[2] * b[0]);
    c[2] = a[0] * b[1] - a[1] * b[0];
    return c;
}

class KDTreeNode {
public:
    KDTreeNode *left;
    KDTreeNode *right;
    int dim;
    real slice;
    std::vector<int> triangles = {};
    bool isLeaf = false;

    KDTreeNode(const pointCloud3d &points, const std::vector<polygonIndexes> &indexesVector,
               std::vector<int>& currTriangles, int depth, int prevSize = -1) {
        // If there is not enough triangles or split is useless, make it a leaf
        if (currTriangles.size() <= 3 || prevSize == currTriangles.size()) {
            isLeaf = true;
            triangles = currTriangles;
            // std::cout << triangles.size() << ", " << depth << std::endl;
            return;
        }

        dim = depth % 3;
        std::vector<real> pnts = {};
        pnts.reserve(currTriangles.size() * 3);
        for (int i: currTriangles) {
            pnts.push_back(points[indexesVector[i][0]][dim]);
            pnts.push_back(points[indexesVector[i][1]][dim]);
            pnts.push_back(points[indexesVector[i][2]][dim]);
        }

        std::sort(pnts.begin(), pnts.end());

        if (pnts.size() % 2 == 0) {
            slice = (pnts[pnts.size() / 2] + pnts[pnts.size() / 2 + 1]) / 2;
        } else {
            slice = pnts[pnts.size() / 2];
        }

        std::vector<int> leftTriangles = {};
        std::vector<int> rightTriangles = {};

        for (int i: currTriangles) {
            bool left = false, right = false;
            for (int j: indexesVector[i]) {
                if (points[j][dim] < slice) {
                    left = true;
                } else {
                    right = true;
                }
            }
            if (left) leftTriangles.push_back(i);
            if (right) rightTriangles.push_back(i);
        }

        left = new KDTreeNode(points, indexesVector, leftTriangles, depth + 1, currTriangles.size());
        right = new KDTreeNode(points, indexesVector, rightTriangles, depth + 1, currTriangles.size());
    }
};

class UnionSet {
    std::vector<int> parent;
    std::vector<int> rank;
    int n_;
public:
    UnionSet(int n) {
        n_ = n;
        parent = {};
        parent.reserve(n);
        rank = {};
        rank.reserve(n);
        makeSet();
    }

    void makeSet() {
        for (int i = 0; i < n_; i++) {
            parent.push_back(i);
            rank.push_back(0);
        }
    }

    int find(int x) {
        if (parent[x] != x) {
            parent[x] = find(parent[x]);
        }

        return parent[x];
    }

    void join(int x, int y) {
        int xset = find(x);
        int yset = find(y);

        if (xset == yset) return;

        if (rank[xset] < rank[yset]) {
            parent[xset] = yset;
        } else if (rank[xset] > rank[yset]) {
            parent[yset] = xset;
        } else {
            parent[yset] = xset;
            rank[xset] = rank[xset] + 1;
        }
    }

    std::vector<int> findAllChildren(int p) {
        std::vector<int> children = {};

        for (int i = 0; i < n_; i++) {
            if (parent[i] == p) children.push_back(i);
        }

        return children;
    }

    std::vector<int> countDistinct() {
        std::vector<int> uniqueParents = {};

        for (int i = 1; i < n_; i++) {
            if (std::find(uniqueParents.begin(), uniqueParents.end(), find(i)) == uniqueParents.end()) {
                uniqueParents.push_back(parent[i]);
            }
        }
        return uniqueParents;
    }

    std::vector<int> countDistinctForSet(std::vector<int> set) {
        std::vector<int> uniqueParents = {};

        for (int i = 1; i < n_; i++) {
            if (std::find(set.begin(), set.end(), i) == set.end()) continue;

            if (std::find(uniqueParents.begin(), uniqueParents.end(), find(i)) == uniqueParents.end()) {
                uniqueParents.push_back(parent[i]);
            }
        }
        return uniqueParents;
    }
};

void getBoundingBox(const pointCloud3d &points, point3d &min, point3d &max) {
    std::vector<real> xpnts = {};
    std::vector<real> ypnts = {};
    std::vector<real> zpnts = {};
    int np = points.size();
    xpnts.reserve(np);
    ypnts.reserve(np);
    zpnts.reserve(np);

    for (const point3d &p: points) {
        xpnts.push_back(p[0]);
        ypnts.push_back(p[1]);
        zpnts.push_back(p[2]);
    }

    std::sort(xpnts.begin(), xpnts.end());
    std::sort(ypnts.begin(), ypnts.end());
    std::sort(zpnts.begin(), zpnts.end());

    min = {xpnts[0], ypnts[0], zpnts[0]};
    max = {xpnts[np - 1], ypnts[np - 1], zpnts[np - 1]};
}

bool rayAABBIntersection(const point3d &orig,
                         const point3d &dir,
                         const point3d &min,
                         const point3d &max,
                         real &tstart,
                         real &tend) {
    real tmin, tmax, t1, t2;

    t1 = (min[0] - orig[0]) / (dir[0] + 1e-8);
    t2 = (max[0] - orig[0]) / (dir[0] + 1e-8);
    if (t1 > t2) {
        real buffor = t1;
        t2 = t1;
        t1 = buffor;
    }

    tmin = t1;
    tmax = t2;

    t1 = (min[1] - orig[1]) / (dir[1] + 1e-8);
    t2 = (max[1] - orig[1]) / (dir[1] + 1e-8);
    if (t1 > t2) {
        real buffor = t1;
        t2 = t1;
        t1 = buffor;
    }

    if (t1 > tmin) tmin = t1;
    if (t2 < tmax) tmax = t2;

    t1 = (min[2] - orig[2]) / (dir[2] + 1e-8);
    t2 = (max[2] - orig[2]) / (dir[2] + 1e-8);
    if (t1 > t2) {
        real buffor = t1;
        t2 = t1;
        t1 = buffor;
    }

    if (t1 > tmin) tmin = t1;
    if (t2 < tmax) tmax = t2;

    if (tmax < 0) {
        return false;
    }

    tstart = 0;
    tend = tmax;
    return true;
}

bool RayTriangleIntersection(const point3d &orig,
                             const point3d &dir,
                             const point3d &v0,
                             const point3d &v1,
                             const point3d &v2,
                             point3d &intersection) {
    // Plane normal
    point3d v0v1 = {v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]};
    point3d v0v2 = {v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]};

    point3d N = crossProduct(v0v1, v0v2);
    real denom = dotProduct(N, N);

    real NdotRayDirection = dotProduct(N, dir);
    if (abs(NdotRayDirection) < 1e-6) {
        return false;
    }

    point3d tmpv0 = {v0[0], v0[1], v0[2]};
    point3d tmporig = {orig[0], orig[1], orig[2]};
    real d = -dotProduct(N, tmpv0);
    real t = -(dotProduct(N, tmporig) + d) / NdotRayDirection;

    if (t < 0) {
        return false;
    }
    point3d P = {orig[0] + t * dir[0], orig[1] + t * dir[1], orig[2] + t * dir[2]};

    point3d C;
    // Edge 0
    point3d edge0 = {v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]};
    point3d vp0 = {P[0] - v0[0], P[1] - v0[1], P[2] - v0[2]};
    C = crossProduct(edge0, vp0);
    if (dotProduct(N, C) < 0) return false;

    // Edge 1
    point3d edge1 = {v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]};
    point3d vp1 = {P[0] - v1[0], P[1] - v1[1], P[2] - v1[2]};
    C = crossProduct(edge1, vp1);
    if (dotProduct(N, C) < 0) return false;

    // Edge 2
    point3d edge2 = {v0[0] - v2[0], v0[1] - v2[1], v0[2] - v2[2]};
    point3d vp2 = {P[0] - v2[0], P[1] - v2[1], P[2] - v2[2]};
    C = crossProduct(edge2, vp2);
    if (dotProduct(N, C) < 0) return false;

    intersection = P;
    return true;
}

struct StackItem {
    KDTreeNode *node;
    real entry, exit;
};

// Traverse kdtree for given ray using recursive traversal algorithm
bool traverseTree(KDTreeNode *rootNode,
                  const point3d &orig,
                  const point3d &dir,
                  const pointCloud3d &points,
                  const std::vector<polygonIndexes> &indexesVector,
                  const point3d &min,
                  const point3d &max,
                  int currTriangleIdx,
                  int &intersectingTriangle) {
    real entry, exit;
    if (!rayAABBIntersection(orig, dir, min, max, entry, exit)) {
        return false;
    }

    std::vector<StackItem> stack = {{rootNode, entry, exit}};

    KDTreeNode *currNode;
    real currEntry, currExit;
    while (stack.size() != 0) {
        StackItem &item = stack[stack.size() - 1];
        currNode = item.node;
        currEntry = item.entry;
        currExit = item.exit;
        stack.pop_back();

        while (!currNode->isLeaf) {
            int a = currNode->dim;
            real t = (currNode->slice - orig[a]) / dir[a];

            KDTreeNode *close;
            KDTreeNode *faraway;
            if (currNode->slice > orig[a]) {
                close = currNode->left;
                faraway = currNode->right;
            } else {
                close = currNode->right;
                faraway = currNode->left;
            }

            if (t > currExit || t < 0) {
                currNode = close;
            } else if (t < currEntry) {
                currNode = faraway;
            } else {
                stack.push_back({ faraway, t, currExit});
                currNode = close;
                currExit = t;
            }
        }

        if (currNode->isLeaf && currNode->triangles.size() > 0) {
            real minDist = 0;
            int closestTriangle = -1;
            point3d p;

            for (int i = 0; i < currNode->triangles.size(); i++) {
                if (currNode->triangles[i] == currTriangleIdx) continue;
                const point3d &v0 = points[indexesVector[currNode->triangles[i]][0]];
                const point3d &v1 = points[indexesVector[currNode->triangles[i]][1]];
                const point3d &v2 = points[indexesVector[currNode->triangles[i]][2]];

                if (RayTriangleIntersection(orig, dir, v0, v1, v2, p)) {
                    real dist = std::sqrt(
                            std::pow(p[0] - orig[0], 2) + std::pow(p[1] - orig[1], 2) + std::pow(p[2] - orig[2], 2));
                    if (closestTriangle < 0 || dist < minDist) {
                        minDist = dist;
                        closestTriangle = currNode->triangles[i];
                    }
                }
            }

            if (closestTriangle > -1) {
                intersectingTriangle = closestTriangle;
                return true;
            }
        }
    }
    return false;
}

real roundReal(real a) {
    return (a * 1e+5);
}

void groupFaces(const char* path, const char* outputPath, int minGroup=1000){
    typedef std::array<real, 3> point3d;
    typedef std::vector<size_t> polygonIndexes;
    typedef std::vector<point3d> pointCloud3d;

    std::string path_(path);
    happly::PLYData input(path_);

    std::vector<polygonIndexes> indexesVector = input.getFaceIndices();
    pointCloud3d points = input.getVertexPositions();

    // Group detected channel faces
    typedef std::tuple<int, int, int> mapKey;
    typedef std::multimap<mapKey, int> multimap;
    typedef multimap::iterator mapIter;

    multimap vertexTriangleConnections;
    
    for(int i=0; i<indexesVector.size(); i++){
        const polygonIndexes& vertices = indexesVector[i];

        // a
        const point3d& pa = points[vertices[0]];
        mapKey keyA = {roundReal(pa[0]), roundReal(pa[1]), roundReal(pa[2])};
        vertexTriangleConnections.insert(std::pair<mapKey, int>(keyA, i));

        // b
        const point3d& pb = points[vertices[1]];
        mapKey keyB = {roundReal(pb[0]), roundReal(pb[1]), roundReal(pb[2])};
        vertexTriangleConnections.insert(std::pair<mapKey, int>(keyB, i));

        //c
        const point3d& pc = points[vertices[2]];
        mapKey keyC = {roundReal(pc[0]), roundReal(pc[1]), roundReal(pc[2])};
        vertexTriangleConnections.insert(std::pair<mapKey, int>(keyC, i));
    }

    mapIter m_it, s_it;
    UnionSet groups(indexesVector.size());

    for (m_it = vertexTriangleConnections.begin();  m_it != vertexTriangleConnections.end(); ){
        mapKey theKey = (*m_it).first;

        std::pair<mapIter, mapIter> keyRange = vertexTriangleConnections.equal_range(theKey);

        for (s_it = keyRange.first;  s_it != keyRange.second;  ++s_it){
            int a = (*m_it).second, b = (*s_it).second;
            if(a == b) continue;
            groups.join(a, b);
        }

        mapIter curr = m_it;
        while (m_it != vertexTriangleConnections.end() && m_it->first == curr->first)
            ++m_it;
    }

    std::vector<int> uniqueParents = groups.countDistinct();
    std::cout << "Detected groups: " << uniqueParents.size() << std::endl;

    // Color each group differently
    // typedef std::array<unsigned char, 3> color;

    // srand(time(NULL));
    // color uniqueColors[uniqueParents.size()] = {};
    // for(int i=0; i<uniqueParents.size(); i++){
    //     uniqueColors[i] = {(unsigned char) (rand() % 256), (unsigned char) (rand() % 256), (unsigned char) (rand() % 256)};
    // }

    // std::vector<std::vector<int>> groupsVertices = {};
    // for(int i=0; i<uniqueParents.size(); i++){
    //     groupsVertices.push_back({});
    // }

    // for(int i=0; i<indexesVector.size(); i++){
    //     auto find = std::find(uniqueParents.begin(), uniqueParents.end(), groups.find(i));
    //     if(find == uniqueParents.end()) continue;
    //     int idx = std::distance(uniqueParents.begin(), find);
    //     groupsVertices[idx].push_back(indexesVector[i][0]);
    //     groupsVertices[idx].push_back(indexesVector[i][1]);
    //     groupsVertices[idx].push_back(indexesVector[i][2]);
    // }

    // std::vector<color> colors = {};

    // for(int i=0; i<points.size(); i++){
    //     int group = -1;
    //     for(int j=0; j<groupsVertices.size(); j++){
    //         if(std::find(groupsVertices[j].begin(), groupsVertices[j].end(), i) != groupsVertices[j].end()){
    //             group = j;
    //             break;
    //         }
    //     }
    //     if(group > -1){
    //         colors.push_back(uniqueColors[group]);
    //     }else{
    //         colors.push_back({255,255,255});
    //     }
    // }

    // input.addVertexColors(colors);
    // input.write("output.ply");    

    std::vector<std::array<double, 3>> meshVertexPositions;
    std::vector<std::vector<size_t>> meshFaceIndices;
    
    for(auto& p : uniqueParents){
        auto group = groups.findAllChildren(p);
        if(group.size() < minGroup) continue;

        for(unsigned long i=0; i<group.size(); i++){
            auto& tri = indexesVector[group[i]];
            meshVertexPositions.push_back(points[tri[2]]);
            meshVertexPositions.push_back(points[tri[1]]);
            meshVertexPositions.push_back(points[tri[0]]);

            std::size_t size = meshVertexPositions.size();

            meshFaceIndices.push_back({size-1, size-2, size-3});
        }
    }

    happly::PLYData outputPLY;

    // Add mesh data (elements are created automatically)
    outputPLY.addVertexPositions(meshVertexPositions);
    outputPLY.addFaceIndices(meshFaceIndices);

    outputPLY.write(outputPath, happly::DataFormat::BinaryBigEndian);
}

int detectInternalChannels(pointCloud3d &points,
                           std::vector<polygonIndexes> &indexesVector,
                           happly::PLYData& input){
    // Create KD Tree
    auto start = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch()
    );
    std::vector<int> currTriangles = {};
    currTriangles.reserve(indexesVector.size());
    for(int i=0; i<indexesVector.size(); i++) currTriangles.push_back(i);

    KDTreeNode* rootNode = new KDTreeNode(points, indexesVector, currTriangles, 0);
    
    auto build = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch()
    );

    // std::cout << "Build tree in: " << build.count() - start.count() << std::endl;

    point3d min, max;
    getBoundingBox(points, min, max);

    // Traverse tree for each face and check if it can be a part of a internal channel
    std::vector<int> possibleChannelFaces = {};
    std::vector<int> intersectingTriangles = {};

    for(int i=0; i<indexesVector.size(); i++){
        const polygonIndexes& face = indexesVector[i];
        const point3d& a = points[face[0]];
        const point3d& b = points[face[1]];
        const point3d& c = points[face[2]];

        point3d normal = crossProduct(b-a, c-a);
        point3d centerPoint = CalculateCenterPointOfTriangle(a, b, c);

        int intersectingTriangle;
        if(traverseTree(rootNode, centerPoint, normal, points, indexesVector, min, max, i, intersectingTriangle)){
            possibleChannelFaces.push_back(i);
            intersectingTriangles.push_back(intersectingTriangle);
        }
    }

    auto traverse = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch()
    );

    // std::cout << "Traversed tree in: " << traverse.count() - build.count() << std::endl;

    // Group detected channel faces
    typedef std::tuple<int, int, int> mapKey;
    typedef std::multimap<mapKey, int> multimap;
    typedef multimap::iterator mapIter;

    multimap vertexTriangleConnections;
    
    for(int i=0; i<possibleChannelFaces.size(); i++){
        const polygonIndexes& vertices = indexesVector[possibleChannelFaces[i]];

        // a
        const point3d& pa = points[vertices[0]];
        mapKey keyA = {roundReal(pa[0]), roundReal(pa[1]), roundReal(pa[2])};
        vertexTriangleConnections.insert(std::pair<mapKey, int>(keyA, i));

        // b
        const point3d& pb = points[vertices[1]];
        mapKey keyB = {roundReal(pb[0]), roundReal(pb[1]), roundReal(pb[2])};
        vertexTriangleConnections.insert(std::pair<mapKey, int>(keyB, i));

        //c
        const point3d& pc = points[vertices[2]];
        mapKey keyC = {roundReal(pc[0]), roundReal(pc[1]), roundReal(pc[2])};
        vertexTriangleConnections.insert(std::pair<mapKey, int>(keyC, i));
    }

    mapIter m_it, s_it;
    UnionSet groups(possibleChannelFaces.size());

    for (m_it = vertexTriangleConnections.begin();  m_it != vertexTriangleConnections.end(); ){
        mapKey theKey = (*m_it).first;

        std::pair<mapIter, mapIter> keyRange = vertexTriangleConnections.equal_range(theKey);

        for (s_it = keyRange.first;  s_it != keyRange.second;  ++s_it){
            int a = (*m_it).second, b = (*s_it).second;
            if(a == b) continue;
            groups.join(a, b);
        }

        mapIter curr = m_it;
        while (m_it != vertexTriangleConnections.end() && m_it->first == curr->first)
            ++m_it;
    }

    // Remove invalid faces
    std::vector<int> channelFaces = {};

    for(int i=0; i<possibleChannelFaces.size(); i++){
        auto t = std::find(possibleChannelFaces.begin(), possibleChannelFaces.end(), intersectingTriangles[i]);
        // If intersecting triangle isn't even a possible channel face, skip
        if(t == possibleChannelFaces.end()) continue;
        // If intersecting triangle is in different group, skip
        if(groups.find(i) != groups.find(std::distance(possibleChannelFaces.begin(), t))) continue;

        channelFaces.push_back(i);
    }

    std::vector<int> uniqueParents = groups.countDistinctForSet(channelFaces);

    auto group = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch()
    );

    // std::cout << "Grouped and removed invalid faces in: " << group.count() - traverse.count() << std::endl;

    // Classify tunels
    real EPSILON = 0.1;
    std::vector<int> invalidTunels = {};

    for(int parent : uniqueParents){
        std::vector<point3d> normals = {};
        for(int channel : channelFaces){
            if(groups.find(channel) == parent){
                const point3d& a = points[indexesVector[possibleChannelFaces[channel]][0]];
                const point3d& b = points[indexesVector[possibleChannelFaces[channel]][1]];
                const point3d& c = points[indexesVector[possibleChannelFaces[channel]][2]];

                point3d normal = crossProduct(b-a, c-a);
                normals.push_back(normal);
            }
        }

        real xBest, yBest, zBest, x, y, z;
        int perpendicularVectors = 0;
        for(int i=100; i>=0; i--){
            real phi = 2*M_PI/i;
            for(int j=0; j<=100; j++){
                real theta = j * M_PI / 100;
                x = std::sin(theta) * std::cos(phi);
                y = std::sin(theta) * std::sin(phi);
                z = std::cos(theta);

                int count = 0;
                for(point3d& normal : normals){
                    real dot = dotProduct(normal, {x,y,z});
                    if(dot < EPSILON && dot > -EPSILON){
                        count++;
                    }
                }

                if(count > perpendicularVectors){
                    perpendicularVectors = count;
                    xBest = x;
                    yBest = y;
                    zBest = z;
                }
            }
        }

        // Check if we can remove tunel from it self by moving in that direction
        bool valid = true;
        for(point3d& normal : normals){
            real dot = dotProduct(normal, {xBest,yBest,zBest});
            if(dot < -EPSILON){
                valid = false;
                break;
            }
        }
        
        if(!valid){
            invalidTunels.push_back(parent);
        }
    }

    auto classify = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch()
    );

    // std::cout << "Classified tunels in: " << classify.count() - group.count() << std::endl;

    // std::cout << "Number of all faces: " << indexesVector.size() << std::endl;
    // std::cout << "Number of channel faces: " << channelFaces.size() << std::endl;
    // std::cout << "Number of detected channels: " << uniqueParents.size() << std::endl;

    // Color each group differently
    // typedef std::array<unsigned char, 3> color;

    // srand(time(NULL));
    // color uniqueColors[invalidTunels.size()] = {};
    // for(int i=0; i<invalidTunels.size(); i++){
    //     uniqueColors[i] = {(unsigned char) (rand() % 256), (unsigned char) (rand() % 256), (unsigned char) (rand() % 256)};
    // }

    // std::vector<std::vector<int>> groupsVertices = {};
    // for(int i=0; i<invalidTunels.size(); i++){
    //     groupsVertices.push_back({});
    // }

    // for(int i=0; i<channelFaces.size(); i++){
    //     auto find = std::find(invalidTunels.begin(), invalidTunels.end(), groups.find(channelFaces[i]));
    //     if(find == invalidTunels.end()) continue;
    //     int idx = std::distance(invalidTunels.begin(), find);
    //     groupsVertices[idx].push_back(indexesVector[possibleChannelFaces[channelFaces[i]]][0]);
    //     groupsVertices[idx].push_back(indexesVector[possibleChannelFaces[channelFaces[i]]][1]);
    //     groupsVertices[idx].push_back(indexesVector[possibleChannelFaces[channelFaces[i]]][2]);
    // }

    // std::vector<color> colors = {};

    // for(int i=0; i<points.size(); i++){
    //     int group = -1;
    //     for(int j=0; j<groupsVertices.size(); j++){
    //         if(std::find(groupsVertices[j].begin(), groupsVertices[j].end(), i) != groupsVertices[j].end()){
    //             group = j;
    //             break;
    //         }
    //     }
    //     if(group > -1){
    //         colors.push_back(uniqueColors[group]);
    //     }else{
    //         colors.push_back({255,255,255});
    //     }
    // }

    // input.addVertexColors(colors);
    // input.write("output.ply");

    auto end = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch()
    );
    // std::cout << "Colored in: " << end.count() - classify.count() << std::endl;

    // std::cout << "Whole process in: " << end.count() - start.count() << std::endl;

    return invalidTunels.size();
}
//*******************************************************************************************
#include <filesystem>
std::string getFileName(std::string path) {
    std::filesystem::path pathToSplit(path);
    std::string fileName = pathToSplit.filename().string();
    int len = fileName.length();
    std::string fileNameWithoutFileFormat = fileName.substr(0, len-4);
    return fileNameWithoutFileFormat;
}

std::string fileNameAfterTransformation(std::string path) {
    std::string fileName = getFileName(path);
    return fileName + "Temp" + ".ply";
}

void writeFile(pointCloud3d& points, std::vector<polygonIndexes>& indexesVector, std::string path) {
    // Create an empty object
    happly::PLYData plyOut;

    // Add mesh data (elements are created automatically)
    plyOut.addVertexPositions(points);
    plyOut.addFaceIndices(indexesVector);
    std::string fileName = fileNameAfterTransformation(path);
    plyOut.write(fileName, happly::DataFormat::ASCII);
}
//*******************************************************************************************

// Polygon Reduction

class Triangle;
class Vertex;

class Triangle {
public:
	Vertex * vertex[3];
	point3d normal;
    Triangle(Vertex *v0, Vertex *v1, Vertex *v2);
	~Triangle();
	void ComputeNormal();
	void ReplaceVertex(Vertex *vold, Vertex *vnew);
	int HasVertex(Vertex *v);
    polygonIndexes getVertices(List<Vertex*>& vertices);
};

class Vertex {
public:
	point3d position;
    int id;
	List<Vertex *> neighbor;
	List<Triangle *> face;
	real objdist;
	Vertex * collapse;
	Vertex(point3d v, int id_);
    ~Vertex();
	void RemoveIfNonNeighbor(Vertex *n);
};

Triangle::Triangle(Vertex *v0,Vertex *v1,Vertex *v2){
	assert(v0!=v1 && v1!=v2 && v2!=v0);  //#mod1
	vertex[0]=v0;
	vertex[1]=v1;
	vertex[2]=v2;
	ComputeNormal();
	for(int i=0;i<3;i++) {
		vertex[i]->face.Add(this);
		for(int j=0;j<3;j++) if(i!=j) {
			vertex[i]->neighbor.AddUnique(vertex[j]);
		}
	}
}

Triangle::~Triangle(){
	int i;
	for(i=0;i<3;i++) {
		if(vertex[i]) vertex[i]->face.Remove(this);
	}
	for(i=0;i<3;i++) {
		int i2 = (i+1)%3;
		if(!vertex[i] || !vertex[i2]) continue;
		vertex[i]->RemoveIfNonNeighbor(vertex[i2]);
		vertex[i2]->RemoveIfNonNeighbor(vertex[i]);
	}
}

int Triangle::HasVertex(Vertex *v) {
	return (v==vertex[0] ||v==vertex[1] || v==vertex[2]);
}

void Triangle::ComputeNormal(){
	point3d v0=vertex[0]->position;
	point3d v1=vertex[1]->position;
	point3d v2=vertex[2]->position;
    point3d a = v1 - v0;
	normal = crossProduct(v1-v0, v2-v1);
	if(magnitude(normal)==0)return;
	normal = normalizeNormal(normal);
}

void Triangle::ReplaceVertex(Vertex *vold,Vertex *vnew) {
	assert(vold && vnew);
	assert(vold==vertex[0] || vold==vertex[1] || vold==vertex[2]);
	assert(vnew!=vertex[0] && vnew!=vertex[1] && vnew!=vertex[2]);
	if(vold==vertex[0]){
		vertex[0]=vnew;
	}
	else if(vold==vertex[1]){
		vertex[1]=vnew;
	}
	else {
		assert(vold==vertex[2]);
		vertex[2]=vnew;
	}
	int i;
	vold->face.Remove(this);
	assert(!vnew->face.Contains(this));
	vnew->face.Add(this);
	for(i=0;i<3;i++) {
		vold->RemoveIfNonNeighbor(vertex[i]);
		vertex[i]->RemoveIfNonNeighbor(vold);
	}
	for(i=0;i<3;i++) {
		assert(vertex[i]->face.Contains(this)==1);
		for(int j=0;j<3;j++) if(i!=j) {
			vertex[i]->neighbor.AddUnique(vertex[j]);
		}
	}
	ComputeNormal();
}

polygonIndexes Triangle::getVertices(List<Vertex*>& vertices){
    return {(unsigned long)vertices.Find(vertex[0]), (unsigned long)vertices.Find(vertex[1]), (unsigned long)vertices.Find(vertex[2])};
}

Vertex::Vertex(point3d v, int id_) {
	position =v;
    id = id_;
}

Vertex::~Vertex(){
	assert(face.num==0);
	while(neighbor.num) {
		neighbor[0]->neighbor.Remove(this);
		neighbor.Remove(neighbor[0]);
	}
	// vertices.Remove(this);
}
void Vertex::RemoveIfNonNeighbor(Vertex *n) {
	if(!neighbor.Contains(n)) return;
	for(int i=0;i<face.num;i++) {
		if(face[i]->HasVertex(n)) return;
	}
	neighbor.Remove(n);
}

real magnitude(point3d v) {
    return (real)sqrt((v[0] * v[0]) + (v[1] * v[1])+ (v[2] * v[2]));
}

point3d normalizeNormal(point3d v) {
    float d=magnitude(v);
    if (d==0) {
		printf("Cant normalize ZERO vector\n");
		assert(0);
		d=0.1f;
	}
    v[0]/=d;
    v[1]/=d;
    v[2]/=d;
    return v;
}

float ComputeEdgeCollapseCost(Vertex *u,Vertex *v) {
	int i;
	float edgelength = magnitude(v->position - u->position);
	float curvature=0;

	List<Triangle *> sides;
	for(i=0;i<u->face.num;i++) {
		if(u->face[i]->HasVertex(v)){
			sides.Add(u->face[i]);
		}
	}
	for(i=0;i<u->face.num;i++) {
		float mincurv=1;
		for(int j=0;j<sides.num;j++) {
			point3d normal1 = u->face[i]->normal;
            point3d normal2 = sides[j]->normal;
            float dotprod = normal1[0]*normal2[0] + normal1[1]*normal2[1] + normal1[2]*normal2[2];
			mincurv = std::min(mincurv,(1-dotprod)/2.0f);
		}
		curvature = std::max(curvature,mincurv);
	}
	return edgelength * curvature;
}

void ComputeEdgeCostAtVertex(Vertex *v) {
	if(v->neighbor.num==0) {
		v->collapse=NULL;
		v->objdist=-0.01f;
		return;
	}
	v->objdist = 1000000;
	v->collapse=NULL;
	for(int i=0;i<v->neighbor.num;i++) {
		float dist;
		dist = ComputeEdgeCollapseCost(v,v->neighbor[i]);
		if(dist<v->objdist) {
			v->collapse=v->neighbor[i];
			v->objdist=dist;
		}
	}
}
void ComputeAllEdgeCollapseCosts(List<Vertex*>& vertices) {
	for(int i=0;i<vertices.num;i++) {
		ComputeEdgeCostAtVertex(vertices[i]);
	}
}

void Collapse(Vertex *u,Vertex *v, List<Vertex*>& vertices, List<Triangle*>& triangles){
    for(int i=0; i<u->face.num; i++){
        if(u->face.num <= i) continue;
        if(u->face[i]->HasVertex(v)){
            triangles.Remove(u->face[i]);
            delete u->face[i];
        }
        u->face[i]->ReplaceVertex(u, v);
    }
    vertices.Remove(u);
    delete u;
}

int getParentVertex(int org, std::map<int, int>& optimizedVertices, std::map<int, int>& pointsVerticesMapper){
    while(optimizedVertices[org] != -1){
        org = optimizedVertices[org];
    }
    return pointsVerticesMapper[org];
}

// Check if collapse is already collapsed to tail of curr
bool isOnTail(int curr, int collapse, std::map<int, int>& optimizedVetrices){
    int i = collapse;
    while(i != -1){
        i = optimizedVetrices[i];
        if(i == curr) return true;
    };
    return false;

}

auto polygonReduction(pointCloud3d& outPoints, 
                      std::vector<polygonIndexes>& outIndexes, 
                      pointCloud3d& inPoints, 
                      std::vector<polygonIndexes>& inIndexes,
                      real COST_EPSILON = 1e-3){ 
    // Join vertices lying on top of each other and create special data structure
    typedef std::tuple<int, int, int> mapKey;
    typedef std::multimap<mapKey, int> multimap;
    typedef multimap::iterator mapIter;
    typedef std::map<int, int> map;

    List<Vertex*> vertices = {};
    List<Triangle*> triangles = {};

    multimap vertexTriangleConnections;
    map deletedVertices;

    auto timer1 = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch()
    );

    for(int i=0; i<inPoints.size(); i++){
        const point3d& p = inPoints[i];
        mapKey key = {roundReal(p[0]), roundReal(p[1]), roundReal(p[2])};
        vertexTriangleConnections.insert(std::pair<mapKey, int>(key, i));
    }

    mapIter m_it, s_it;

    for (m_it = vertexTriangleConnections.begin();  m_it != vertexTriangleConnections.end(); ){
        mapKey theKey = (*m_it).first;

        std::pair<mapIter, mapIter> keyRange = vertexTriangleConnections.equal_range(theKey);

        Vertex *v = new Vertex(inPoints[(*m_it).second], vertices.num);
        vertices.Add(v);

        for (s_it = keyRange.first;  s_it != keyRange.second;  ++s_it){
            int b = (*s_it).second;
            deletedVertices[b] = vertices.num - 1;
        }

        mapIter curr = m_it;
        while (m_it != vertexTriangleConnections.end() && m_it->first == curr->first)
            ++m_it;
    }

    for(int i=0; i<inIndexes.size(); i++){
        int a = deletedVertices[inIndexes[i][0]];
        int b = deletedVertices[inIndexes[i][1]];
        int c = deletedVertices[inIndexes[i][2]];

        Triangle *t = new Triangle(vertices[a], vertices[b], vertices[c]);
        triangles.Add(t);
    }
    
    ComputeAllEdgeCollapseCosts(vertices);

    std::map<int, int> optimizedVertices;
    std::map<int, int>::iterator it;
    std::map<int, int> pointsVerticesMapper;

    // Collapsed cannot be already collapsed to the tail of i !!
    for(int i=0; i<vertices.num; i++){
        if(vertices[i]->objdist <= COST_EPSILON && !isOnTail(vertices[i]->id, vertices[i]->collapse->id, optimizedVertices)){
            optimizedVertices[vertices[i]->id] = vertices[i]->collapse->id;
        }else{
            optimizedVertices[vertices[i]->id] = -1;
        }
    }

    for(it = optimizedVertices.begin(); it != optimizedVertices.end(); it++){
        if(it->second == -1){
            outPoints.push_back(vertices[it->first]->position);
            pointsVerticesMapper[vertices[it->first]->id] = outPoints.size() - 1;
        }
    }

    for(int i=0; i<triangles.num; i++){
        unsigned long v0 = getParentVertex(triangles[i]->vertex[0]->id, optimizedVertices, pointsVerticesMapper);
        unsigned long v1 = getParentVertex(triangles[i]->vertex[1]->id, optimizedVertices, pointsVerticesMapper);
        unsigned long v2 = getParentVertex(triangles[i]->vertex[2]->id, optimizedVertices, pointsVerticesMapper);
        // Triangle has been deleted
        if(v0 == v1 || v0 == v2 || v1 == v2) continue;

        outIndexes.push_back({v0, v1, v2});
    }

    auto timer2 = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch()
    );

    return timer2.count() - timer1.count();
}

static real distanceBetweenPointAndSegment(point3d& bg, point3d& end, point3d& p){
    point3d ab = {end[0] - bg[0], end[1] - bg[1], end[2] - bg[2]};
    point3d av = {p[0] - bg[0], p[1] - bg[1], p[2] - bg[2]};

    if(dotProduct(av, ab) <= 0){
        return sqrt(av[0]*av[0] + av[1]*av[1] + av[2]*av[2]);
    }

    point3d bv = {p[0] - end[0], p[1] - end[1], p[2] - end[2]};

    if(dotProduct(bv, ab) >= 0){
        return sqrt(bv[0]*bv[0] + bv[1]*bv[1] + bv[2]*bv[2]);
    }

    point3d cross = crossProduct(ab, av);
    return sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]) / sqrt(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
}

// Remove connectors using markers glued to every connector
void removeConnectorsWithMarkers(const char* path, const char* outputPath, std::vector<std::pair<point3d, point3d>>& markers){
    typedef std::array<real, 3> point3d;
    typedef std::vector<size_t> polygonIndexes;
    typedef std::vector<point3d> pointCloud3d;
    typedef std::tuple<int, int, int> mapKey;
    typedef std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, int>>> sharedMap;
    typedef sharedMap::iterator mapIter;

    auto timer1 = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch()
    );

    std::string path_(path);
    happly::PLYData input(path_);

    std::vector<polygonIndexes> indexesVector = input.getFaceIndices();
    pointCloud3d points = input.getVertexPositions();

    std::vector<polygonIndexes> sharedIndexesVector = {};
    pointCloud3d sharedPoints = {};

    sharedMap uniqueVertices;

    real markerPos = 0.5;

    const real THRESH = 10;

    for(int i=0; i<points.size(); i++){
        int a = roundReal(points[i][0]); 
        int b = roundReal(points[i][1]); 
        int c = roundReal(points[i][2]);

        if(uniqueVertices[a][b][c] == 0){
            uniqueVertices[a][b][c] = sharedPoints.size() + 1;
            sharedPoints.push_back(points[i]);
        }
    }

    auto timer2 = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch()
    );

    std::cout << timer2.count() - timer1.count() << std::endl;

    pointCloud3d withoutConnectors = {};
    std::map<int, int> map;

    timer1 = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch()
    );

    for(auto& pair : markers){
        point3d vec = {pair.second[0] - pair.first[0], pair.second[1] - pair.first[1], pair.second[2] - pair.first[2]};
        
        point3d dist = {vec[0] * markerPos, vec[1] * markerPos, vec[2] * markerPos};
        point3d bg = {pair.first[0] - dist[0], pair.first[1] - dist[1], pair.first[2] - dist[2]};
        point3d end = {pair.second[0] + dist[0], pair.second[1] + dist[1], pair.second[2] + dist[2]};

        for(int i=0; i<sharedPoints.size(); i++){
            real distanceFromConnector = distanceBetweenPointAndSegment(bg, end, sharedPoints[i]);

            if(distanceFromConnector > THRESH){
                withoutConnectors.push_back(sharedPoints[i]);
                map[i+1] = withoutConnectors.size() - 1;
            }else{
                map[i+1] = -1;
            }
        }
    }

    for(auto& tri : indexesVector){
        unsigned long a = uniqueVertices[roundReal(points[tri[0]][0])][roundReal(points[tri[0]][1])][roundReal(points[tri[0]][2])];
        unsigned long b = uniqueVertices[roundReal(points[tri[1]][0])][roundReal(points[tri[1]][1])][roundReal(points[tri[1]][2])];
        unsigned long c = uniqueVertices[roundReal(points[tri[2]][0])][roundReal(points[tri[2]][1])][roundReal(points[tri[2]][2])];

        if(map[a] == -1 || map[b] == -1 || map[c] == -1) continue;

        sharedIndexesVector.push_back({(unsigned long)map[a],(unsigned long)map[b],(unsigned long)map[c]});
    }

    timer2 = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch()
    );

    std::cout << timer2.count() - timer1.count() << std::endl;

    happly::PLYData outputPLY;

    // Add mesh data (elements are created automatically)
    outputPLY.addVertexPositions(withoutConnectors);
    outputPLY.addFaceIndices(sharedIndexesVector);

    outputPLY.write(outputPath, happly::DataFormat::BinaryBigEndian);
}

auto normalizeVector(point3d vector_3d) {
    if (vector_3d[0] == 0 && vector_3d[1] == 0 && vector_3d[2] == 0) {
        return point3d({ 0, 0, 0 });
    }
    real vector_length = sqrt(pow(vector_3d[0], 2) + pow(vector_3d[1], 2) + pow(vector_3d[2], 2));
    for (int i = 0; i < 3; i++) {
        vector_3d[i] = vector_3d[i] / vector_length;
    }
    return vector_3d;
}

// Remove connectors knowing only position of the frame, corners are top left, top right, bottom right, bottom left
void removeConnectorsWithCorners(const char* path, const char* outputPath, std::array<point3d, 4>& corners){
    typedef std::array<real, 3> point3d;
    typedef std::vector<size_t> polygonIndexes;
    typedef std::vector<point3d> pointCloud3d;
    typedef std::tuple<int, int, int> mapKey;
    typedef std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, int>>> sharedMap;
    typedef sharedMap::iterator mapIter;

    // We must very carefully adjust connector length because it can remove part of model (not sure how to do it yet).
    const double artificalConnectorLen = 100;

    // Threshold that controls how close vertex should be to line segment to be classified as artificial connector
    const double DIST_THRESH = 10;

    // For how many positions on frame edge algorithm will be run
    const double N_EDGE_SEGMENTS = 20;

    for(int i=0; i<4; i++){
        std::cout << i << ", " << (i+1)%4 << std::endl;
        // Calculate line segment between corners (it will start in i'th corner)
        point3d mainLineVec = corners[(i+1)%4] - corners[i];

        double stepSize = std::sqrt(mainLineVec[0]*mainLineVec[0] + mainLineVec[1]*mainLineVec[1] + mainLineVec[2]*mainLineVec[2]) / N_EDGE_SEGMENTS;

        point3d mainLineVecNorm = normalizeVector(mainLineVec);
        // Calculate directional vector of artifical connector
        point3d artificialConnectorVec = normalizeVector(corners[(i+2)%4] - corners[(i+1)%4]);
        artificialConnectorVec[0] *= artificalConnectorLen;
        artificialConnectorVec[1] *= artificalConnectorLen;
        artificialConnectorVec[2] *= artificalConnectorLen;

        // Maximum number of faces around
        int max = 0;

        for(int j=0; j<N_EDGE_SEGMENTS; j++){
            point3d artificialConnectorStart = corners[i];
            artificialConnectorStart[0] += mainLineVecNorm[0] * (j+stepSize);
            artificialConnectorStart[1] += mainLineVecNorm[1] * (j+stepSize);
            artificialConnectorStart[2] += mainLineVecNorm[2] * (j+stepSize);

            // Calculate how many vertices are around line segment and pick the one with most vertices
        }
    }

}