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
    std::vector<ui> max_degrees;
    std::vector<std::vector<VertexID>> vert;
    std::vector<std::vector<ui>> bin;
    std::vector<std::vector<ui>> pos;
    std::vector<std::vector<ui>> degrees;
    std::vector<std::vector<std::vector<VertexID>>> adj;
public:
    XYCore();
    XYCore& operator=(const XYCore& other);
    void xyCoreInitialization(Graph& graph, bool sort = false);
    void generateXYCore(const Graph &graph, Graph &subgraph, ui x, ui y, bool is_exact, bool is_map, bool is_copy);

    ui getDelta(const Graph &graph);
    ui skyline_core_num(Graph &graph, ui cur, ui x, ui y, bool reduced = false);
private:
    inline void decDeg(ui cur, VertexID t, std::vector<std::vector<VertexID>> &vert_copy, std::vector<std::vector<ui>> &bin_copy, std::vector<std::vector<ui>> &pos_copy,
                       std::vector<std::vector<ui>> &degrees_copy);
};


#endif //DENSESTSUBGRAPH_XYCORE_H
