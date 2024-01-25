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
    std::vector<std::vector<VertexID>> vert;
    std::vector<std::vector<ui>> bin;
    std::vector<std::vector<ui>> pos;
    std::vector<std::vector<ui>> degrees;
public:
    XYCore();
    void xyCoreInitialization(const Graph &graph);
    void generateXYCore(const Graph &graph, Graph &x_y_core, ui x, ui y, bool is_exact = true);
    ui getDelta(const Graph &graph);
    ui skyline_core_num(Graph &graph, ui cur, ui x, ui y, bool reduced = false);
private:
    inline void decDeg(ui cur, VertexID t, std::vector<std::vector<VertexID>> &vert_copy, std::vector<std::vector<ui>> &bin_copy, std::vector<std::vector<ui>> &pos_copy,
                       std::vector<std::vector<ui>> &degrees_copy);
};


#endif //DENSESTSUBGRAPH_XYCORE_H
