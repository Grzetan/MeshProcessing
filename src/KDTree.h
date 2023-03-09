#pragma once
#include <vector>
#include <algorithm>
#include "types.h"

class KDTreeNode {
public:
    KDTreeNode *left;
    KDTreeNode *right;
    int dim;
    real slice;
    std::vector<int> triangles = {};
    bool isLeaf = false;

    KDTreeNode(pointCloud3d &points, const std::vector<polygonIndexes> &indexesVector,
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