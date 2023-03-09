#include <iostream>
#include <math.h>

struct Vec3{
    double x, y, z;

    Vec3 operator+(Vec3& v){
        return {x + v.x, y + v.y, z + v.z};
    }

    Vec3 operator-(Vec3& v){
        return {x - v.x, y - v.y, z - v.z};
    }

    Vec3 operator*(double v){
        return {x * v, y * v, z * v};
    }

    Vec3 operator/(double v){
        return {x / v, y / v, z / v};
    }

    double& operator[](unsigned int idx){
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

    Vec3 crossProduct(Vec3& v){
        return {
            y * v.z - z * v.y,
            x * v.z - z * v.x,
            x * v.y - y * v.x
        };
    }

    double dotProduct(Vec3& v){
        return x * v.x + y * v.y + z * v.z;
    }

    double magnitude(){
        return std::sqrt(x*x + y*y + z*z);
    }

    double distance(Vec3& v){
        return std::sqrt(std::pow(x - v.x, 2) + std::pow(y - v.y, 2) + std::pow(z - v.z, 2));
    }

    friend std::ostream& operator<<(std::ostream& os, const Vec3& dt){
        return os << dt.x << ", " << dt.y << ", " << dt.z << std::endl;
    }
};

int main(){
    // Right Down = 61.77, -142.32, −561.81
    // Right Up = 14.82, 139.63, −711.93
    // Left Down = -10.35, -37.94, −335.5
    // Left Up = -60.98, 247.93, −482.51

    Vec3 rightDown{74.6, -114.42, -569.60};
    Vec3 rightUp{32.54, 123.32, -695.90};
    Vec3 upLeft{-38.6, 241.1, -501.9};

    Vec3 topVector = upLeft - rightUp;

    Vec3 downLeft = rightDown + topVector;

    std::cout << downLeft; // -> 3.46, 3.36, -375.6
}