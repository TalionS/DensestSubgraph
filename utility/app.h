#ifndef DENSESTSUBGRAPH_APP_H
#define DENSESTSUBGRAPH_APP_H

#include <vector>
#include <iostream>
#include "types.h"
#include "graph.h"
#include <fstream>
#include <string>
#include <algorithm>
#include <random>
#include <cmath>



class CoreApp{

public:
    ui nodes_count;
    ui pos;
    ui k;
    ui size;
    double opt;
    std::vector<ui> w;
    std::vector<bool> selected;
    std::vector<ui> id;
    std::vector<ui> *new_edge;
    std::vector<ui> deg;


public:

    explicit CoreApp(ui nodes_count = 0, ui pos = 0);

    void Init(Graph &graph);

};

class PKMC{

public:
    ui hmax;
    ui s;
    std::vector<ui> h;
    explicit PKMC();
    void Init(Graph &graph);
};

class Greedy{

public:
    std::vector<double> h;
    void Init(Graph &graph);
};
#endif //DENSESTSUBGRAPH_FLOWNETWORK_H