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
    std::vector<ui> w;


public:

    explicit CoreApp(ui nodes_count = 0, ui pos = 0);

    void Init(Graph &graph);

};

#endif //DENSESTSUBGRAPH_FLOWNETWORK_H