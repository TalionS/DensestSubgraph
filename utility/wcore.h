//
// Created by Talion on 2024/1/11.
//

#include <vector>
#include <algorithm>
#include <climits>
#include "types.h"
#include "graph.h"

#ifndef DENSESTSUBGRAPH_WCORE_H
#define DENSESTSUBGRAPH_WCORE_H

class WCore{
public:
    ui w;
    std::vector<std::pair<VertexID, VertexID>> edges;
    std::vector<ui> induce_numbers;
    std::vector<VertexID> vertices[2];
    std::vector<ui> degrees[2];
public:
    void generateMaxWCore(Graph &graph, Graph &subgraph);
    void wCoreDecomposition(Graph &graph);
    void getMaxCNPair(Graph &graph, std::pair<ui, ui> &max_core_num_pair);
};

#endif //DENSESTSUBGRAPH_WCORE_H
