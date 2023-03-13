#include "types.h"
#include "functions.h"

int main(){
    std::string path = "./objects/1merged.ply";

    std::array<point3d, 4> corners;
    corners[0] = {-60.98, 247.93, -482.51}; // Top left
    corners[1] = {14.82, 139.63, -711.93}; // Top right
    corners[2] = {61.77, -142.32, -561.81}; // Bottom right
    corners[3] = {-10.35, -37.94, -335.5}; // Bottom left

    removeConnectorsWithCorners(path.c_str(), "testowanie.ply", corners);

    // std::vector<std::pair<point3d, point3d>> markers = { { {74.6, -114.42, -569.60}, {73.6, -113.42, -568.60} }};
    // std::vector<std::pair<point3d, point3d>> markers = { { {14.02, -20.3, -486}, {24, -80, -454} }};

    // removeConnectorsWithMarkers(path.c_str(), "testowanie.ply", markers);
}