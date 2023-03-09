#pragma once 
#include <chrono>
#include "types.h"
#include "KDTree.h"
#include "UnionSet.h"
#include "happly.h"

point3d CalculateCenterPointOfTriangle(point3d &a, point3d &b, point3d &c) {
    point3d center = {(a[0] + b[0] + c[0]) / 3,
                      (a[1] + b[1] + c[1]) / 3,
                      (a[2] + b[2] + c[2]) / 3};
    return center;
}

void getBoundingBox(pointCloud3d &points, point3d &min, point3d &max) {
    std::vector<real> xpnts = {};
    std::vector<real> ypnts = {};
    std::vector<real> zpnts = {};
    int np = points.size();
    xpnts.reserve(np);
    ypnts.reserve(np);
    zpnts.reserve(np);

    for (point3d &p: points) {
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

bool rayAABBIntersection(point3d &orig,
                         point3d &dir,
                         point3d &min,
                         point3d &max,
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

bool RayTriangleIntersection(point3d &orig,
                             point3d &dir,
                             point3d &v0,
                             point3d &v1,
                             point3d &v2,
                             point3d &intersection) {
    // Plane normal
    point3d v0v1 = {v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]};
    point3d v0v2 = {v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]};

    point3d N = v0v1.crossProduct(v0v2);
    real denom = N.dotProduct(N);

    real NdotRayDirection = N.dotProduct(dir);
    if (abs(NdotRayDirection) < 1e-6) {
        return false;
    }

    point3d tmpv0 = {v0[0], v0[1], v0[2]};
    point3d tmporig = {orig[0], orig[1], orig[2]};
    real d = -N.dotProduct(tmpv0);
    real t = -(N.dotProduct(tmporig) + d) / NdotRayDirection;

    if (t < 0) {
        return false;
    }
    point3d P = {orig[0] + t * dir[0], orig[1] + t * dir[1], orig[2] + t * dir[2]};

    point3d C;
    // Edge 0
    point3d edge0 = {v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]};
    point3d vp0 = {P[0] - v0[0], P[1] - v0[1], P[2] - v0[2]};
    C = edge0.crossProduct(vp0);
    if (N.dotProduct(C) < 0) return false;

    // Edge 1
    point3d edge1 = {v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]};
    point3d vp1 = {P[0] - v1[0], P[1] - v1[1], P[2] - v1[2]};
    C = edge1.crossProduct(vp1);
    if (N.dotProduct(C) < 0) return false;

    // Edge 2
    point3d edge2 = {v0[0] - v2[0], v0[1] - v2[1], v0[2] - v2[2]};
    point3d vp2 = {P[0] - v2[0], P[1] - v2[1], P[2] - v2[2]};
    C = edge2.crossProduct(vp2);
    if (N.dotProduct(C) < 0) return false;

    intersection = P;
    return true;
}

struct StackItem {
    KDTreeNode *node;
    real entry, exit;
};

// Traverse kdtree for given ray using recursive traversal algorithm
bool traverseTree(KDTreeNode *rootNode,
                  point3d &orig,
                  point3d &dir,
                  pointCloud3d &points,
                  const std::vector<polygonIndexes> &indexesVector,
                  point3d &min,
                  point3d &max,
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
                point3d &v0 = points[indexesVector[currNode->triangles[i]][0]];
                point3d &v1 = points[indexesVector[currNode->triangles[i]][1]];
                point3d &v2 = points[indexesVector[currNode->triangles[i]][2]];

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
        point3d& pa = points[vertices[0]];
        mapKey keyA = {roundReal(pa[0]), roundReal(pa[1]), roundReal(pa[2])};
        vertexTriangleConnections.insert(std::pair<mapKey, int>(keyA, i));

        // b
        point3d& pb = points[vertices[1]];
        mapKey keyB = {roundReal(pb[0]), roundReal(pb[1]), roundReal(pb[2])};
        vertexTriangleConnections.insert(std::pair<mapKey, int>(keyB, i));

        //c
        point3d& pc = points[vertices[2]];
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

    std::vector<point3d> meshVertexPositions;
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
        point3d& a = points[face[0]];
        point3d& b = points[face[1]];
        point3d& c = points[face[2]];

        point3d normal = (b-a).crossProduct(c-a);
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
        point3d& pa = points[vertices[0]];
        mapKey keyA = {roundReal(pa[0]), roundReal(pa[1]), roundReal(pa[2])};
        vertexTriangleConnections.insert(std::pair<mapKey, int>(keyA, i));

        // b
        point3d& pb = points[vertices[1]];
        mapKey keyB = {roundReal(pb[0]), roundReal(pb[1]), roundReal(pb[2])};
        vertexTriangleConnections.insert(std::pair<mapKey, int>(keyB, i));

        //c
        point3d& pc = points[vertices[2]];
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
                point3d& a = points[indexesVector[possibleChannelFaces[channel]][0]];
                point3d& b = points[indexesVector[possibleChannelFaces[channel]][1]];
                point3d& c = points[indexesVector[possibleChannelFaces[channel]][2]];

                point3d normal = (b-a).crossProduct(c-a);
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
                    real dot = normal.dotProduct({x,y,z});
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
            real dot = normal.dotProduct({xBest,yBest,zBest});
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

real distanceBetweenPointAndSegment(point3d& bg, point3d& end, point3d& p){
    point3d ab = {end[0] - bg[0], end[1] - bg[1], end[2] - bg[2]};
    point3d av = {p[0] - bg[0], p[1] - bg[1], p[2] - bg[2]};

    if(av.dotProduct(ab) <= 0){
        return sqrt(av[0]*av[0] + av[1]*av[1] + av[2]*av[2]);
    }

    point3d bv = {p[0] - end[0], p[1] - end[1], p[2] - end[2]};

    if(bv.dotProduct(ab) >= 0){
        return sqrt(bv[0]*bv[0] + bv[1]*bv[1] + bv[2]*bv[2]);
    }

    point3d cross = ab.crossProduct(av);
    return sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]) / sqrt(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
}

// Remove connectors using markers glued to every connector
void removeConnectorsWithMarkers(const char* path, const char* outputPath, std::vector<std::pair<point3d, point3d>>& markers){
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

    pointCloud3d withoutConnectors = {};
    std::map<int, int> map;

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