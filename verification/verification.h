//
// Created by yy on 12/1/23.
//

#ifndef DENSESTSUBGRAPH_VERIFICATION_H
#define DENSESTSUBGRAPH_VERIFICATION_H

#include "../utility/types.h"
#include "vector"
#include "../utility/graph.h"
#include "../utility/flownetwork.h"
#include "../utility/lp.h"
#include "../utility/xycore.h"
#include <cmath>
#include <algorithm>

class Verification{
public:
    bool flowExactVerification(Graph &graph, XYCore xy_core, double l, double r);
    bool directedLPExactVerification(Graph &graph, Graph &x_y_core, LinearProgramming &lp, std::pair<ui, ui> best_pos, std::vector<std::vector<VertexID>> &vertices, std::pair<double, double> ratios, double rho_c);
    bool UndirectedflowExactVerification(Graph &graph, double l, double r);
    bool UndirectedlpVerification(Graph &graph, LinearProgramming &lp, FlowNetwork &flow, std::vector<VertexID> *vertices);
};

#endif //DENSESTSUBGRAPH_VERIFICATION_H
