#include "types.h"
#include "functions.h"

int main(){
    std::string path = "./objects/1merged.ply";

    std::vector<std::pair<point3d, point3d>> markers = { { {74.6, -114.42, -569.60}, {73.6, -113.42, -568.60} }};
    removeConnectorsWithMarkers(path.c_str(), "testowanie.ply", markers);
}