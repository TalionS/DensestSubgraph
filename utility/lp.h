#ifndef DENSESTSUBGRAPH_LP_H
#define DENSESTSUBGRAPH_LP_H

#include <vector>
#include <iostream>
#include "types.h"
#include "graph.h"
#include <fstream>
#include <string>
#include <algorithm>
#include <random>
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

public:
    ui nodes_count_;
    ui edges_count_;
    ui type_;
    ui sort_type;
    ui cur_iter_num;
public:
    std::vector<std::vector<double>> r;
    std::vector<Alpha> alpha;
    std::vector<Alpha> beta;
    std::vector<double> weight;


public:
    ~LinearProgramming();
    explicit LinearProgramming(bool is_directed, ui type = 0, ui vertices_count = 0, ui edge_count = 0, ui sort = 0);

    void Iterate(double learning_rate, double ratio = 0, bool is_synchronous = false);

    void FistaIterate(double learning_rate, double t, double ratio = 0, bool is_synchronous = false);

    void MWUIterate(ui t, bool is_synchronous = false);

    void Init(Graph &graph, double ratio = 0);

    void sort(Graph &graph);
};

#endif //DENSESTSUBGRAPH_FLOWNETWORK_H