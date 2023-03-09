#pragma once
#include <vector>
#include <algorithm>

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