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
void Interface::initQuery(float *qx, float *qy, float *qz, size_t chainidx1 ) {
    chainIdx1 = chainidx1;
    for(unsigned int i = 0; i < queryLength; i++) {
        query_coordinates[i][0] = qx[i];
        query_coordinates[i][1] = qy[i];
        query_coordinates[i][2] = qz[i];
    }
    query_grid = Grid(query_coordinates, queryLength, true);
}

void Interface::getinterface(unsigned int targetLen, float *tx, float *ty, float *tz, std::map<unsigned int, std::vector<unsigned int>> &qInterfaceIndex, size_t chainidx2 ) {
    unsigned int targetLength = targetLen;
    unsigned int chainIdx2 = chainidx2;
    float **target_coordinates = new float*[targetLength];
    for(unsigned int i = 0; i < targetLength; i++) {
        target_coordinates[i] = new float[3];
    }
    for(unsigned int i = 0; i < targetLength; i++) {
        target_coordinates[i][0] = tx[i];
        target_coordinates[i][1] = ty[i];
        target_coordinates[i][2] = tz[i];
    }
    Interface::Grid target_grid = Grid(target_coordinates, targetLength, false);
    const int DIR = 14;
    int dx[DIR] = {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0};
    int dy[DIR] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1, 1, 1, 1, 0};
    int dz[DIR] = {0, 1, 0, -1, 1, 0, -1, 1, 0, -1, 1, 0, -1, 1};
    std::vector<std::tuple<int, int, int>> tboxes;
    std::map<std::tuple<int, int, int>, bool> visited_boxes;
    for (unsigned int target_idx = 0; target_idx < targetLength; target_idx++){
        std::tuple<int, int, int> tbox_coord = target_grid.getGridCoordinates(target_coordinates[target_idx]);
        tboxes.push_back(tbox_coord);
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
                        std::vector<unsigned int>& indexes1 = qInterfaceIndex[chainIdx1];
                        if (std::find(indexes1.begin(), indexes1.end(), index1) == indexes1.end()) {
                            qInterfaceIndex[chainIdx1].push_back(index1);
                        }
                    }
                    for (size_t i2 = neighbor_members.first; i2 < neighbor_members.second; i2++) {
                        index2 = target_grid.box[i2].second;
                        std::vector<unsigned int>& indexes2 = qInterfaceIndex[chainIdx2];
                        if (std::find(indexes2.begin(), indexes2.end(), index2) == indexes2.end()) {
                            qInterfaceIndex[chainIdx2].push_back(index2);
                        }
                    }
                }
            }
        }
    }
    if(target_coordinates) {
        for(unsigned int i = 0; i < targetLength; i++) {
            delete[] target_coordinates[i];
        }
        delete[] target_coordinates;
    }
    tboxes.clear();
    visited_boxes.clear();
}