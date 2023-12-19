//
// Created by yy on 12/16/23.
//

#include <vector>
#include <algorithm>
#include "types.h"
#include "graph.h"

#ifndef DENSESTSUBGRAPH_XYCORE_H
#define DENSESTSUBGRAPH_XYCORE_H

class XYCore{
public:
    std::vector<VertexID> vert[2];
    std::vector<ui> bin[2];
    std::vector<ui> pos[2];
    std::vector<ui> degrees[2];
public:
    XYCore();
    void xyCoreInitialization(const Graph &graph);
    void generateXYCore(const Graph &graph, Graph &x_y_core, ui x, ui y);
    ui getDelta(const Graph &graph);
    ui skyline_core_num(Graph &graph, ui cur, ui x, ui y, bool reduced = false);
private:
    inline void decDeg(ui cur, VertexID t);
};


#endif //DENSESTSUBGRAPH_XYCORE_H
