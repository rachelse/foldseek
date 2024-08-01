#include "MultimerInterface.h"
#include <string.h>
#include <algorithm>

Interface::Interface(float cutoff) :cutoff(cutoff) {}

Interface::~Interface() {
    if(query_coordinates) {
        for(unsigned int i = 0; i < queryLength; i++) {
            delete[] query_coordinates[i];
        }
        delete[] query_coordinates;
    }
    if(target_coordinates) {
        for(unsigned int i = 0; i < targetLength; i++) {
            delete[] target_coordinates[i];
        }
        delete[] target_coordinates;
    }
}
void Interface::initQuery(unsigned int queryLength, float *qx, float *qy, float *qz, size_t chainidx1 ) {
    queryLength = queryLength;
    chainIdx1 = chainidx1;
    query_coordinates = new float*[queryLength];
    for(unsigned int i = 0; i < queryLength; i++) {
        query_coordinates[i] = new float[3];
    }
    for(unsigned int i = 0; i < queryLength; i++) {
        query_coordinates[i][0] = qx[i];
        query_coordinates[i][1] = qy[i];
        query_coordinates[i][2] = qz[i];
    }
    query_grid = Grid(query_coordinates, queryLength);
}

void Interface::getinterface(unsigned int targetLen, float *tx, float *ty, float *tz, std::vector<std::set<unsigned int>> &qInterfaceIndex, size_t chainidx2 ) {
    targetLength = targetLen;
    chainIdx2 = chainidx2;
    target_coordinates = new float*[targetLength];
    for(unsigned int i = 0; i < targetLength; i++) {
        target_coordinates[i] = new float[3];
    }
    for(unsigned int i = 0; i < targetLength; i++) {
        target_coordinates[i][0] = tx[i];
        target_coordinates[i][1] = ty[i];
        target_coordinates[i][2] = tz[i];
    }
    target_grid = Grid(target_coordinates, targetLength);
    const int DIR = 14;
    int dx[DIR] = {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0};
    int dy[DIR] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1, 1, 1, 1, 0};
    int dz[DIR] = {0, 1, 0, -1, 1, 0, -1, 1, 0, -1, 1, 0, -1, 1};
    std::vector<std::tuple<int, int, int>> tboxes;
    std::map<std::tuple<int, int, int>, bool> visited_boxes;
    for (unsigned int target_idx = 0; target_idx < targetLength){
        std::tuple<int, int, int> tbox_coord = target_grid.getGridCoordinates(target_coordinates[target_idx]);
        tboxes.push_back(tbox_coord);
    }
    for (unsigned int query_idx = 0; query_idx < queryLength){
        std::tuple<int, int, int> qbox_coord = query_grid.getGridCoordinates(query_coordinates[query_idx]);
        if (visited_boxes.find(qbox_coord) != visited_boxes.end()) {
            continue;
        }
        visited_boxes[qbox_coord] = true;
        std::pair<size_t, size_t> box_members = query_grid.getBoxMemberRange(box_coord);
        for (int dir = 0; dir < DIR; dir++) {
            std::tuple<int, int, int> key = std::make_tuple(std::get<0>(box_coord) + dx[dir], std::get<1>(box_coord) + dy[dir], std::get<2>(box_coord) + dz[dir]);
            std::pair<size_t, size_t> neighbor_members = target_grid.getBoxMemberRange(key);
            if (neighbor_members.second != 0){
                for (size_t i = box_members.first; i < box_members.second; i++){
                    int index1 = query_grid.box[i].second;
                    if (qInterfaceIndex[chainIdx1].find(index1) == qInterfaceIndex[chainIdx1].end()){
                        qInterfaceIndex[chainIdx1].insert(index1);
                    }
                }
                size_t i2 = (dir == 0) ? i + 1 : neighbor_members.first;
                for (; i2 < neighbor_members.second; i2++) {
                    int index2 = target_grid.box[i2].second;
                    if (qInterfaceIndex[chainIdx2].find(index2) == qInterfaceIndex[chainIdx2].end()){
                        qInterfaceIndex[chainIdx2].insert(index2);
                    }
                }
            }
        }
    }
}