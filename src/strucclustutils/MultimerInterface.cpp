#include "MultimerInterface.h"
#include <string.h>
#include <algorithm>

Interface::Interface(unsigned int maxQueryLength, unsigned int maxTargetLength, float cutoff)
        : maxQueryLength(maxQueryLength), maxTargetLength(maxTargetLength) {
    maxAlignLength = std::max(maxQueryLength, maxTargetLength);
    query_coordinates = new float*[maxQueryLength];
    for(unsigned int i = 0; i < maxQueryLength; i++) {
        query_coordinates[i] = new float[3];
    }
    target_coordinates = new float*[maxTargetLength];
    for(unsigned int i = 0; i < maxTargetLength; i++) {
        target_coordinates[i] = new float[3];
    }
    dists_to_score = new bool*[maxQueryLength];
    for(unsigned int i = 0; i < maxQueryLength; i++) {
        dists_to_score[i] = new bool[maxQueryLength];
    }
    query_to_align = new int[maxQueryLength];
    target_to_align = new int[maxTargetLength];
    align_to_query = new int[maxAlignLength];
    align_to_target = new int[maxAlignLength];
}

Interface::~Interface() {
    if(query_coordinates) {
        for(unsigned int i = 0; i < maxQueryLength; i++) {
            delete[] query_coordinates[i];
        }
        delete[] query_coordinates;
    }
    if(target_coordinates) {
        for(unsigned int i = 0; i < maxTargetLength; i++) {
            delete[] target_coordinates[i];
        }
        delete[] target_coordinates;
    }
    if(dists_to_score) {
        for(unsigned int i = 0; i < maxQueryLength; i++) {
            delete[] dists_to_score[i];
        }
        delete[] dists_to_score;
    }
}

void Interface::initQuery(unsigned int queryLen, float *qx, float *qy, float *qz) {
    queryLength = queryLen;
    for(unsigned int i = 0; i < queryLength; i++) {
        query_coordinates[i][0] = qx[i];
        query_coordinates[i][1] = qy[i];
        query_coordinates[i][2] = qz[i];
    }
    // Initialize arrays
    for(unsigned int i = 0; i < queryLength; i++) {
        memset(dists_to_score[i], 0, sizeof(bool) * queryLength);
    }

    query_grid = Grid(query_coordinates, queryLength);

    for(unsigned int col = 0; col < queryLength; col++) {
        for (unsigned int row = 0; row < queryLength; row++) {
            float distance = BasicFunction::dist(query_coordinates[row], query_coordinates[col]);
            bool isClose = (col != row) && (distance < cutoff*cutoff);
            dists_to_score[col][row] = isClose;
            dists_to_score[row][col] = isClose;
        }
    }

}

void Interface::getinterface(unsigned int targetLen, float *tx, float *ty, float *tz, std::vector<std::set<unsigned int>> &qInterfaceIndex) {
    targetLength = targetLen;

    for(unsigned int i = 0; i < targetLength; i++) {
        target_coordinates[i][0] = tx[i];
        target_coordinates[i][1] = ty[i];
        target_coordinates[i][2] = tz[i];
    }
    //sooyoung 0731 TODO
    //get interface region using DIR, box like below
    //return: vector<set> of indexes of interface..(like qInterfacePos) 
}

// void Interface::calculateLddtScores() {
//     const int DIR = 14;
//     int dx[DIR] = {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0};
//     int dy[DIR] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1, 1, 1, 1, 0};
//     int dz[DIR] = {0, 1, 0, -1, 1, 0, -1, 1, 0, -1, 1, 0, -1, 1};
//     std::map<std::tuple<int, int, int>, bool> visited_boxes;
//     // Create a set of unique boxes containing aligned residues
//     for (unsigned int align_idx = 0; align_idx < alignLength; align_idx++) {
//         if (align_to_query[align_idx] != -1) {
//             std::tuple<int, int, int> box_coord = query_grid.getGridCoordinates(query_coordinates[align_to_query[align_idx]]);
//             // If box_coord has been visited already, skip
//             if (visited_boxes.find(box_coord) != visited_boxes.end()) {
//                 continue;
//             }
//             visited_boxes[box_coord] = true;
//             std::pair<std::tuple<int, int, int>, int> ref;
//             ref.first = box_coord;
//             ref.second = 0;

//             std::pair<size_t, size_t> box_members = query_grid.getBoxMemberRange(box_coord);

//             for (size_t i = box_members.first; i < box_members.second; i++) {
//                 int query_idx1 = query_grid.box[i].second;
//                 int align_idx1 = query_to_align[query_idx1];
//                 if (align_idx1 == -1) {
//                     continue;
//                 }

//                 // Different boxes
//                 for (int dir = 0; dir < DIR; dir++) {
//                     std::pair<std::tuple<int, int, int>, int> ref;
//                     std::tuple<int, int, int> key = std::make_tuple(std::get<0>(box_coord) + dx[dir],
//                                     std::get<1>(box_coord) + dy[dir],
//                                     std::get<2>(box_coord) + dz[dir]);
//                     std::pair<size_t, size_t> boxPrime_members = query_grid.getBoxMemberRange(key);

//                     size_t i2 = (dir == 0) ? i + 1 : boxPrime_members.first;
//                     for (; i2 < boxPrime_members.second; i2++) {
//                         int query_idx2 = query_grid.box[i2].second;;
//                         int align_idx2 = query_to_align[query_idx2];
//                         if (align_idx2 == -1) {
//                             continue;
//                         }
//                         if (dists_to_score[query_idx1][query_idx2]) {
//                             float distance = BasicFunction::dist(query_coordinates[query_idx1], query_coordinates[query_idx2]);
//                             float dist_sub = BasicFunction::dist(target_coordinates[align_to_target[align_idx1]],
//                                                   target_coordinates[align_to_target[align_idx2]]);
//                             float d_l = std::abs(distance - dist_sub);
//                             // float score = 0.25 * ((d_l < 0.5*0.5) + (d_l < 1.0*1.0) + (d_l < 2.0*2.0) + (d_l < 4.0*4.0));
//                             // reduce_score[align_idx2] += score;
//                             // reduce_score[align_idx1] += score;
//                         }
//                     }
//                 }
//             }
//         }
//     }
// }
