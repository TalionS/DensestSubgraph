//
// Created by yy on 12/1/23.
//

#ifndef DENSESTSUBGRAPH_EXTRACTION_H
#define DENSESTSUBGRAPH_EXTRACTION_H

#include "../utility/types.h"
#include "vector"
#include "../utility/graph.h"
#include "../utility/flownetwork.h"
#include "../utility/lp.h"
#include <cmath>
#include <algorithm>
#include "wcore.h"

class Extraction{
public:
    void flowExactExtraction(Graph &graph, FlowNetwork &flow, double &l, double &r, ui &s_size, ui &t_size);
    void directedBSApproExtraction(Graph &graph, std::vector<std::vector<bool>> is_peeled, std::vector<std::vector<VertexID>> &vertices);
    void directedPMApproExtraction(Graph &graph, ui edges_count, std::vector<std::vector<VertexID>> vertices);
    void directedCPExtraction(Graph &graph, LinearProgramming &lp, std::pair<ui, ui> &best_pos,
                              std::vector<std::vector<VertexID>> &vertices, std::pair<double, double> ratios,
                              double &ratio_o, double &ratio_p, double &rho, double &rho_c);
    void directedVWApproExtraction(Graph &graph, LinearProgramming &lp,
                                   std::vector<std::vector<VertexID>> &vertices, std::pair<double, double> ratios,
                                   double &rho, double &vw_rho);
    void directedWCoreApproExtraction(Graph &graph, WCore &w_core);
    void UndirectedflowExactExtraction(Graph &graph, FlowNetwork &flow, double &l, double &r, std::vector<VertexID> *vertices);
    void UndirectedlpExactExtraction(Graph &graph, LinearProgramming &lp, std::vector<VertexID> *vertices);
};

#endif //DENSESTSUBGRAPH_EXTRACTION_H
