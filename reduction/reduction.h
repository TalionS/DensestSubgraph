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
#include <climits>
#include "lp.h"
#include "wcore.h"

class Reduction {
public:
    void xyCoreReduction(Graph &graph, Graph &x_y_core, std::pair<double, double> ratios, double &l, double &r,
                         bool &is_init, bool is_dc, bool is_divide_by_number, bool is_exact, bool is_map,
                         bool is_res, ui res_width, bool is_copy);
    void kCoreReduction(Graph &graph, double &l, double &r);
    void stableSetReduction(Graph &graph, LinearProgramming &lp, std::vector<std::pair<VertexID, VertexID>> &edges, bool stable_set_reduction, bool is_map = false);
    void wCoreReduction(Graph &graph, Graph &subgraph, WCore &w_core);
    void UndirectedkCoreReduction(Graph &graph, LinearProgramming &lp, bool weighted = false);
    void UndirectedStableReduction(Graph &graph, LinearProgramming &lp);

public:
    std::vector<ui> core;
    XYCore xycore;


private:
    void coreDecomposition(const Graph &graph);

};

#endif //DENSESTSUBGRAPH_REDUCTION_H
