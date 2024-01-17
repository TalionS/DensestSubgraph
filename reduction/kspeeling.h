//
// Created by Talion on 2024/1/15.
//

#include "types.h"
#include "graph.h"
#include <vector>
#include <boost/heap/fibonacci_heap.hpp>
#ifndef DENSESTSUBGRAPH_KSPEELING_H
#define DENSESTSUBGRAPH_KSPEELING_H
class KSPeeling {
private:
    std::vector<std::vector<VertexID>> vertices_;
    std::vector<std::vector<VertexID>> vert_;
    std::vector<std::vector<ui>> bin_;
    std::vector<std::vector<ui>> pos_;
    std::vector<std::vector<ui>> degrees_;
private:
    void peeling(double &density, std::vector<std::vector<VertexID>> &vertices);
    inline void decDeg(ui cur, VertexID t, std::vector<std::vector<VertexID>> &vert_copy, std::vector<std::vector<ui>> &bin_copy, std::vector<std::vector<ui>> &pos_copy,
                       std::vector<std::vector<ui>> &degrees_copy);
public:
    KSPeeling(Graph &graph);
};

#endif //DENSESTSUBGRAPH_KSPEELING_H
