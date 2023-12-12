//
// Created by yy on 12/1/23.
//

#ifndef DENSESTSUBGRAPH_VERIFICATION_H
#define DENSESTSUBGRAPH_VERIFICATION_H

#include "../utility/types.h"
#include "vector"
#include "../utility/graph.h"
#include "../utility/flownetwork.h"
#include <cmath>

class Verification{
public:
    bool flowExactVerification(Graph &graph, double l, double r);
    bool UndirectedflowExactVerification(Graph &graph, double l, double r);
};

#endif //DENSESTSUBGRAPH_VERIFICATION_H
