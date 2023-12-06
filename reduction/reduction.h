//
// Created by yy on 23-11-15.
//

#ifndef DENSESTSUBGRAPH_REDUCTION_H
#define DENSESTSUBGRAPH_REDUCTION_H

#include "utility/graph.h"
#include "utility/types.h"

#include <queue>
#include <vector>
#include <iostream>
#include <list>
#include <algorithm>
#include <cmath>

class Reduction {
public:
    void xyCoreReduction(Graph &graph, Graph &x_y_core, std::pair<double, double> ratio, double &l, double &r, bool &is_init, bool is_dc);

public:
    std::vector<ui> core;

    std::vector<VertexID> *vert;
    std::vector<ui> *pos;
    std::vector<ui> *bin;
private:
    void coreDecomposition(const Graph &graph);
    void generateXYCore(const Graph &graph, Graph &x_y_core, ui x, ui y);
    void xyCoreInitialization(const Graph &graph);

};

#endif //DENSESTSUBGRAPH_REDUCTION_H
