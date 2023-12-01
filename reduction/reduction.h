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

class Reduction {
public:
    std::vector<ui> coreDecomposition(const Graph &graph);

    std::vector<VertexID> *xyCoreDecomposition(const Graph &graph, ui x, ui y);

};

#endif //DENSESTSUBGRAPH_REDUCTION_H
