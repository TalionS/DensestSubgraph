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
#include "../utility/app.h"
#include <cmath>
#include <algorithm>

class Verification{
public:
    bool flowExactVerification(Graph &graph, double l, double r);
    bool directedBSApproVerification(Graph &graph, ui edges_count, std::vector<std::vector<VertexID>> vertices);
    bool directedFixedKSApproVerification(Graph &graph, ui &cnt, ui edges_count,
                                          std::vector<std::vector<VertexID>> vertices);
    bool directedCPVerification(Graph &graph, Graph &subgraph, LinearProgramming &lp, std::pair<ui, ui> best_pos,
                                std::vector<std::vector<VertexID>> &vertices, std::pair<double, double> ratios,
                                double rho,
                                double rho_c, double &ratio_o, double &ratio_p, bool &stable_set_reduction,
                                std::vector<std::pair<VertexID, VertexID>> &edges, double epsilon = 0);
    bool directedPMApproVerification(std::vector<std::vector<VertexID>> vertices);
    bool directedVWApproVerification(Graph &graph,
                                     LinearProgramming &lp,
                                     std::vector<std::vector<VertexID>> &vertices,
                                     double rho,
                                     double vw_rho,
                                     double epsilon = 0);
    bool UndirectedflowExactVerification(Graph &graph, double l, double r);
    bool UndirectedlpVerification(Graph &graph, LinearProgramming &lp, FlowNetwork &flow, std::vector<VertexID> *vertices);
    bool UndirectedCoreAppVerification(Graph &graph, CoreApp &ca);
};

#endif //DENSESTSUBGRAPH_VERIFICATION_H
