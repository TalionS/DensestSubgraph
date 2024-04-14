#ifndef DENSESTSUBGRAPH_EXTRACTION_H
#define DENSESTSUBGRAPH_EXTRACTION_H

#include "../utility/types.h"
#include "vector"
#include "../utility/graph.h"
#include "../utility/flownetwork.h"
#include "../utility/lp.h"
#include "xycore.h"
#include <cmath>
#include <algorithm>
#include "wcore.h"

class Extraction{
public:
    void flowExactExtraction(Graph &graph, Graph &subgraph, std::pair<double, double> ratio, FlowNetwork &flow,
                             double &l, double &r, double &ratio_o, double &ratio_p, bool is_map);
    void directedCoreApproExtraction(Graph &graph, Graph &subgraph, std::pair<ui, ui> &max_core_num_pair);
    void directedBSApproExtraction(Graph &graph, std::vector<std::vector<bool>> is_peeled, std::vector<std::vector<VertexID>> &vertices);
    void directedPMApproExtraction(Graph &graph, ui edges_count, std::vector<std::vector<VertexID>> vertices);
    void directedCPExtraction(Graph &graph, LinearProgramming &lp, std::pair<ui, ui> &best_pos,
                              std::vector<std::vector<VertexID>> &vertices, std::pair<double, double> ratios,
                              double &ratio_o, double &ratio_p, double &rho, double &rho_c, bool is_map);
//    void directedVWApproExtraction(Graph &graph, LinearProgramming &lp,
//                                   std::vector<std::vector<VertexID>> &vertices, std::pair<double, double> ratios,
//                                   double &rho, double &vw_rho);
    void directedVWApproExtraction(Graph &graph, Graph &vw_graph, std::vector<std::vector<VertexID>> vertices,
                                   std::vector<VertexID> *verticess, double &ratio_o, double &ratio_p);
    void UndirectedflowExactExtraction(Graph &graph, FlowNetwork &flow, double &l, double &r, std::vector<VertexID> *vertices);
    void UndirectedlpExactExtraction(Graph &graph, LinearProgramming &lp, std::vector<VertexID> *vertices);
};

#endif //DENSESTSUBGRAPH_EXTRACTION_H
