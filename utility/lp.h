#ifndef DENSESTSUBGRAPH_LP_H
#define DENSESTSUBGRAPH_LP_H

#include <vector>
#include <iostream>
#include "types.h"
#include "graph.h"
#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>


struct Alpha{
    ui id_first;
    ui id_second;
    double weight_first;
    double weight_second;
};

class LinearProgramming{

private:
    bool is_directed_;
    ui nodes_count_;
    ui edges_count_;
public:
    std::vector<double> *r;
    std::vector<Alpha> *alpha;

public:
    ~LinearProgramming();
    explicit LinearProgramming(bool is_directed, ui vertices_count = 0, ui edge_count = 0);

    void Iterate(double learning_rate, double ratio = 0);

    void Init(Graph &graph, double ratio = 0);

};

#endif //DENSESTSUBGRAPH_FLOWNETWORK_H