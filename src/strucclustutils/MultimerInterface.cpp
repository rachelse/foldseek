#include "MultimerInterface.h"
#include <string.h>
#include <algorithm>

float Interface::min[3] = { Interface::INF, Interface::INF, Interface::INF };

Interface::Interface(unsigned int queryLength) :queryLength(queryLength) {
    query_coordinates1 = new float*[queryLength];
    target_coordinates1 = new float*[queryLength];
    for(unsigned int i = 0; i < queryLength; i++) {
        query_coordinates1[i] = new float[3];
        target_coordinates1[i] = new float[3];
    }
}

Interface::~Interface() {
    if(query_coordinates1) {
        for(unsigned int i = 0; i < queryLength; i++) {
            delete[] query_coordinates1[i];
            delete[] target_coordinates1[i];
        }
        delete[] query_coordinates1;
        delete[] target_coordinates1;
    }
}
void Interface::initQuery(float *qx1, float *qy1, float *qz1, float *tx1, float *ty1, float *tz1, unsigned int chainidx1 ) {
    for(unsigned int i = 0; i < queryLength; i++) {
        query_coordinates1[i][0] = qx1[i];
        query_coordinates1[i][1] = qy1[i];
        query_coordinates1[i][2] = qz1[i];
        target_coordinates1[i][0] = tx1[i];
        target_coordinates1[i][1] = ty1[i];
        target_coordinates1[i][2] = tz1[i];
    }
    query_grid1 = Grid(query_coordinates1, queryLength, true);
    chainidx = chainidx1;
}

void Interface::getinterface(unsigned int targetLen, float *qx2, float *qy2, float *qz2, float *tx2, float *ty2, float *tz2, AlignedCoordinate &qInterface, AlignedCoordinate &tInterface, unsigned int chainidx2) {
    unsigned int targetLength = targetLen;
    float **query_coordinates2 = new float*[targetLength];
    float **target_coordinates2 = new float*[targetLength];
    for(unsigned int i = 0; i < targetLength; i++) {
        query_coordinates2[i] = new float[3];
        target_coordinates2[i] = new float[3];
    }
    for(unsigned int i = 0; i < targetLength; i++) {
        query_coordinates2[i][0] = qx2[i];
        query_coordinates2[i][1] = qy2[i];
        query_coordinates2[i][2] = qz2[i];
        target_coordinates2[i][0] = tx2[i];
        target_coordinates2[i][1] = ty2[i];
        target_coordinates2[i][2] = tz2[i];
    }
    Interface::Grid query_grid2;
    query_grid2 = Grid(query_coordinates2, targetLength, false);
    const int DIR = 19;
    int dx[DIR] = {0, 0, 0, 0, 0, 0, 0,  0,  0,  1,  1, -1, -1, 1, 1, -1, -1, 1, -1};
    int dy[DIR] = {0, 0, 1, 0, -1, 1, 1, -1, -1, -1, 1, -1, 1, 0, 0, 0, 0, 0, 0};
    int dz[DIR] = {0, 1, 0, -1, 0, -1, 1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0};
    std::vector<std::tuple<int, int, int>> tboxes;
    std::map<std::tuple<int, int, int>, bool> visited_boxes;
    for (unsigned int target_idx = 0; target_idx < targetLength; target_idx++){
        std::tuple<int, int, int> tbox_coord = query_grid2.getGridCoordinates(query_coordinates2[target_idx]);
        if (std::find(tboxes.begin(), tboxes.end(), tbox_coord) == tboxes.end()){
            tboxes.push_back(tbox_coord);
        }
    }
    for (unsigned int query_idx = 0; query_idx < queryLength; query_idx++){
        std::tuple<int, int, int> qbox_coord = query_grid1.getGridCoordinates(query_coordinates1[query_idx]);
        if (visited_boxes.find(qbox_coord) != visited_boxes.end()) {
            continue;
        }
        visited_boxes.emplace(qbox_coord, true);
        std::pair<size_t, size_t> box_members = query_grid1.getBoxMemberRange(qbox_coord);
        for (int dir = 0; dir < DIR; dir++) {
            std::tuple<int, int, int> key = std::make_tuple(std::get<0>(qbox_coord) + dx[dir], std::get<1>(qbox_coord) + dy[dir], std::get<2>(qbox_coord) + dz[dir]);
            if (std::find(tboxes.begin(), tboxes.end(), key) != tboxes.end()){
                std::pair<size_t, size_t> neighbor_members = query_grid2.getBoxMemberRange(key);
                unsigned int index1, index2;
                if (neighbor_members.second-neighbor_members.first > 0){
                    for (size_t i = box_members.first; i < box_members.second; i++){
                        index1 = query_grid1.box[i].second;
                        std::vector<unsigned int>& indexes1 = qInterface.chainResidueIndexMap[chainidx];
                        if (std::find(indexes1.begin(), indexes1.end(), index1) == indexes1.end()) {
                            qInterface.chainResidueIndexMap[chainidx].push_back(index1);
                            qInterface.x.push_back(query_coordinates1[index1][0]);
                            qInterface.y.push_back(query_coordinates1[index1][1]);
                            qInterface.z.push_back(query_coordinates1[index1][2]);
                            tInterface.x.push_back(target_coordinates1[index1][0]);
                            tInterface.y.push_back(target_coordinates1[index1][1]);
                            tInterface.z.push_back(target_coordinates1[index1][2]);
                        }
                    }
                    for (size_t i2 = neighbor_members.first; i2 < neighbor_members.second; i2++) {
                        index2 = query_grid2.box[i2].second;
                        std::vector<unsigned int>& indexes2 = qInterface.chainResidueIndexMap[chainidx2];
                        if (std::find(indexes2.begin(), indexes2.end(), index2) == indexes2.end()) {
                            qInterface.chainResidueIndexMap[chainidx2].push_back(index2);
                            qInterface.x.push_back(query_coordinates2[index2][0]);
                            qInterface.y.push_back(query_coordinates2[index2][1]);
                            qInterface.z.push_back(query_coordinates2[index2][2]);
                            tInterface.x.push_back(target_coordinates2[index2][0]);
                            tInterface.y.push_back(target_coordinates2[index2][1]);
                            tInterface.z.push_back(target_coordinates2[index2][2]);
                        }
                    }
                }
            }
        }
    }
    for(unsigned int i = 0; i < targetLength; i++) {
        delete[] query_coordinates2[i];
        delete[] target_coordinates2[i];
    }
    delete[] query_coordinates2;
    delete[] target_coordinates2;
    tboxes.clear();
    visited_boxes.clear();
}