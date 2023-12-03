//
// Created by yy on 12/1/23.
//

#ifndef DENSESTSUBGRAPH_ALLOCATION_H
#define DENSESTSUBGRAPH_ALLOCATION_H

#include "../utility/flownetwork.h"
#include "../utility/graph.h"

class Allocation{
public:
    void flowExactAllocation(Graph &graph, FlowNetwork &flow, double ratio, double *l, double *r);
};

#endif //DENSESTSUBGRAPH_ALLOCATION_H
