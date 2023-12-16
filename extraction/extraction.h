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

class Extraction{
public:
    void flowExactExtraction(Graph &graph, FlowNetwork &flow, double &l, double &r, std::vector<VertexID> *vertices);
    void UndirectedflowExactExtraction(Graph &graph, FlowNetwork &flow, double &l, double &r, std::vector<VertexID> *vertices);
    void UndirectedlpExactExtraction(Graph &graph, LinearProgamming &lp, std::vector<VertexID> *vertices);
};

#endif //DENSESTSUBGRAPH_EXTRACTION_H
