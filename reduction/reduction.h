//
// Created by yy on 23-11-15.
//

#ifndef DENSESTSUBGRAPH_REDUCTION_H
#define DENSESTSUBGRAPH_REDUCTION_H

#include "utility/graph.h"
#include "utility/types.h"
#include "utility/xycore.h"

#include <queue>
#include <vector>
#include <iostream>
#include <list>
#include <algorithm>
#include <cmath>

class Reduction {
public:
    void xyCoreReduction(Graph &graph, Graph &x_y_core, std::pair<double, double> ratio, double &l, double &r, bool &is_init, bool is_dc);
    void kCoreReduction(Graph &graph, double &l, double &r);

public:
    std::vector<ui> core;
    XYCore xycore;


private:
    void coreDecomposition(const Graph &graph);

};

#endif //DENSESTSUBGRAPH_REDUCTION_H
