#include "MultimerInterface.h"
#include <string.h>
#include <algorithm>

float Interface::min[3] = { Interface::INF, Interface::INF, Interface::INF };

Interface::Interface(unsigned int queryLength) :queryLength(queryLength) {
    query_coordinates = new float*[queryLength];
    for(unsigned int i = 0; i < queryLength; i++) {
        query_coordinates[i] = new float[3];
    }
}

Interface::~Interface() {
    if(query_coordinates) {
        for(unsigned int i = 0; i < queryLength; i++) {
            delete[] query_coordinates[i];
        }
        delete[] query_coordinates;
    }
}
void Interface::initQuery(float *qx, float *qy, float *qz ) {
    for(unsigned int i = 0; i < queryLength; i++) {
        query_coordinates[i][0] = qx[i];
        query_coordinates[i][1] = qy[i];
        query_coordinates[i][2] = qz[i];
    }
    query_grid = Grid(query_coordinates, queryLength, true);
}

void Interface::getinterface(unsigned int targetLen, float *tx, float *ty, float *tz, AlignedCoordinate &qInterface) {
    unsigned int targetLength = targetLen;
    float **target_coordinates = new float*[targetLength];
    for(unsigned int i = 0; i < targetLength; i++) {
        target_coordinates[i] = new float[3];
    }
    for(unsigned int i = 0; i < targetLength; i++) {
        target_coordinates[i][0] = tx[i];
        target_coordinates[i][1] = ty[i];
        target_coordinates[i][2] = tz[i];
    }
    Interface::Grid target_grid;
    target_grid = Grid(target_coordinates, targetLength, false);
    const int DIR = 19;
    int dx[DIR] = {0, 0, 0, 0, 0, 0, 0,  0,  0,  1,  1, -1, -1, 1, 1, -1, -1, 1, -1};
    int dy[DIR] = {0, 0, 1, 0, -1, 1, 1, -1, -1, -1, 1, -1, 1, 0, 0, 0, 0, 0, 0};
    int dz[DIR] = {0, 1, 0, -1, 0, -1, 1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0};
    std::vector<std::tuple<int, int, int>> tboxes;
    std::map<std::tuple<int, int, int>, bool> visited_boxes;
    for (unsigned int target_idx = 0; target_idx < targetLength; target_idx++){
        std::tuple<int, int, int> tbox_coord = target_grid.getGridCoordinates(target_coordinates[target_idx]);
        if (std::find(tboxes.begin(), tboxes.end(), tbox_coord) == tboxes.end()){
            tboxes.push_back(tbox_coord);
        }
    }
    for (unsigned int query_idx = 0; query_idx < queryLength; query_idx++){
        std::tuple<int, int, int> qbox_coord = query_grid.getGridCoordinates(query_coordinates[query_idx]);
        if (visited_boxes.find(qbox_coord) != visited_boxes.end()) {
            continue;
        }
        visited_boxes.emplace(qbox_coord, true);
        std::pair<size_t, size_t> box_members = query_grid.getBoxMemberRange(qbox_coord);
        for (int dir = 0; dir < DIR; dir++) {
            std::tuple<int, int, int> key = std::make_tuple(std::get<0>(qbox_coord) + dx[dir], std::get<1>(qbox_coord) + dy[dir], std::get<2>(qbox_coord) + dz[dir]);
            if (std::find(tboxes.begin(), tboxes.end(), key) != tboxes.end()){
                std::pair<size_t, size_t> neighbor_members = target_grid.getBoxMemberRange(key);
                unsigned int index1, index2;
                if (neighbor_members.second-neighbor_members.first > 0){
                    for (size_t i = box_members.first; i < box_members.second; i++){
                        index1 = query_grid.box[i].second;
                        qInterface.x.push_back(query_coordinates[index1][0]);
                        qInterface.y.push_back(query_coordinates[index1][1]);
                        qInterface.z.push_back(query_coordinates[index1][2]);
                    }
                    for (size_t i2 = neighbor_members.first; i2 < neighbor_members.second; i2++) {
                        index2 = target_grid.box[i2].second;
                        qInterface.x.push_back(target_coordinates[index2][0]);
                        qInterface.y.push_back(target_coordinates[index2][1]);
                        qInterface.z.push_back(target_coordinates[index2][2]);
                    }
                }
            }
        }
    }
    for(unsigned int i = 0; i < targetLength; i++) {
        delete[] target_coordinates[i];
    }
    delete[] target_coordinates;
    tboxes.clear();
    visited_boxes.clear();
}