// #include <iostream>
#include <string>
#include <unordered_map>
#include <cmath>
#include <map>
#include <limits>
#include <vector>
#include "FastSort.h"
#include "Debug.h"

#ifndef MULTIMERINTERFACE_H
#define MULTIMERINTERFACE_H

class Interface {
public:
    static constexpr float CUTOFF = 8.0;
    static constexpr float INF = std::numeric_limits<float>::infinity();
    static float min[3];

    struct Grid {
        Grid() {};
        Grid(float **& m1, unsigned int queryLength, bool isInit) {
            int len = queryLength;
            if (isInit == true){
                for(int i = 0; i < len; i++) {
                    for(int dim = 0; dim < 3; dim++) {
                        if(m1[i][dim] < Interface::min[dim]) Interface::min[dim] = m1[i][dim];
                    }
                }
            }
            box.clear();
            box_index.clear();
            for(int i = 0; i < (int)queryLength; i++) {
                int box_coord[3];
                for(int dim = 0; dim < 3; dim++) {
                    box_coord[dim] = (int)((m1[i][dim] - Interface::min[dim]) / CUTOFF);                    
                }
                box.emplace_back(std::make_tuple(box_coord[0], box_coord[1], box_coord[2]), i);
            }
            
            SORT_SERIAL(box.begin(), box.end());
            // precompute index for box to get range per entry
            // tuple -> (start, end) in box_index
            std::tuple<int, int, int> currBoxKey = box[0].first;
            size_t currBoxStart = 0;
            for(size_t i = 0; i < box.size(); i++){
                if(box[i].first != currBoxKey){
                    box_index.push_back(std::make_pair(currBoxKey, std::make_pair(currBoxStart, i)));
                    currBoxKey = box[i].first;
                    currBoxStart = i;
                }
            }
            box_index.push_back(std::make_pair(currBoxKey, std::make_pair(currBoxStart, box.size())));
        }

        std::tuple<int, int, int> getGridCoordinates(const float *point) {
            int box_coord[3];
            for(int dim = 0; dim < 3; dim++) {
                box_coord[dim] = (int)((point[dim] - Interface::min[dim]) / CUTOFF);
            }
            return std::make_tuple(box_coord[0], box_coord[1], box_coord[2]);
        }

        struct PairTupleComparator {
            bool operator()(const std::pair<std::tuple<int, int, int>, std::pair<long unsigned int, long unsigned int>>& lhs,
                            const std::pair<std::tuple<int, int, int>, std::pair<long unsigned int, long unsigned int>>& rhs) const {
                // Compare the tuples first
                if (lhs.first < rhs.first) {
                    return true;
                }
                if (rhs.first < lhs.first) {
                    return false;
                }

                // If the tuples are equal, compare the pairs
                return lhs.second < rhs.second;
            }
        };

        std::pair<size_t, size_t> getBoxMemberRange(std::tuple<int, int, int> &box_coord) {
            auto it = std::lower_bound(box_index.begin(), box_index.end(), std::make_pair(box_coord, std::make_pair(0, 0)), Grid::PairTupleComparator());
            if(it == box_index.end() || it->first != box_coord) {
                return std::make_pair(0, 0);
            }
            return it->second;
        }
        std::vector<std::pair<std::tuple<int, int, int>, int>> box;
        std::vector<std::pair<std::tuple<int, int, int>, std::pair<size_t, size_t>>> box_index;
    };
    Interface(unsigned int queryLen);
    ~Interface();

    void initQuery(float *qx, float *qy, float *qz, size_t chainidx1 );
    void getinterface(unsigned int targetLen, float *tx, float *ty, float *tz, std::map<unsigned int, std::vector<unsigned int>> &qInterfaceIndex, size_t chainidx2);

private:
    unsigned int chainIdx1;
    unsigned int queryLength;
    float **query_coordinates;
    Interface::Grid query_grid;
};

#endif
