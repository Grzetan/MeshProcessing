#pragma once
#include <iostream>
#include <array>
#include <vector>
#include <map>
#include <unordered_map>
#include <math.h>

typedef double real;

struct point3d{
    real x, y, z;

    point3d operator+(point3d v){
        return {x + v.x, y + v.y, z + v.z};
    }

    point3d operator+=(point3d v){
        return {x + v.x, y + v.y, z + v.z};
    }

    point3d operator-(point3d v){
        return {x - v.x, y - v.y, z - v.z};
    }

    point3d operator-=(point3d v){
        return {x - v.x, y - v.y, z - v.z};
    }

    point3d operator*(double v){
        return {x * v, y * v, z * v};
    }

    point3d operator/(double v){
        return {x / v, y / v, z / v};
    }

    real& operator[](int idx){
        if(idx < 0 || idx > 2) throw std::runtime_error("Index out of range");
        switch(idx){
            case 0:
                return x;
                break;
            case 1:
                return y;
                break;
            case 2:
                return z;
                break;
            default:
                return x;
                break;
        }
    }

    point3d crossProduct(point3d v){
        return {
            y * v.z - z * v.y,
            x * v.z - z * v.x,
            x * v.y - y * v.x
        };
    }

    double dotProduct(point3d v){
        return x * v.x + y * v.y + z * v.z;
    }

    double magnitude(){
        return std::sqrt(x*x + y*y + z*z);
    }

    double distance(point3d& v){
        return std::sqrt(std::pow(x - v.x, 2) + std::pow(y - v.y, 2) + std::pow(z - v.z, 2));
    }

    friend std::ostream& operator<<(std::ostream& os, const point3d& dt){
        return os << "{ " << dt.x << ", " << dt.y << ", " << dt.z << " }";
    }
};

typedef std::vector<size_t> polygonIndexes;
typedef std::vector<point3d> pointCloud3d;
typedef std::tuple<int, int, int> mapKey;
typedef std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, int>>> sharedMap;
typedef sharedMap::iterator mapIter;