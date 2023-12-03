//
// Created by yy on 12/1/23.
//

#include "extraction.h"

void Extraction::flowExactExtraction(Graph &graph, FlowNetwork &flow, double *l, double *r, std::vector<VertexID> *vertices) {
    std::vector<VertexID> tmp_S;
    ui n = graph.getVerticesCount();
    VertexID s = 0, t = 2 * n + 1;
    vertices[0].clear();
    vertices[1].clear();
    double mid = (*l + *r) / 2.0;

    flow.getMinCut(s, t, tmp_S);
    ui edge_num = 0;
    for (auto &v: tmp_S) {
        if (v == 0) continue;
        if (v <= n) {
            vertices[0].push_back(v - 1);
            for (auto &edge: flow.adj_[v]) {
                if (edge.to > n && edge.to <= 2 * n && flow.dist_[edge.to] >= n) {
                    edge_num++;
                }
            }
        } else vertices[1].push_back(v - n - 1);
    }

    if (!vertices[0].empty() && !vertices[1].empty()) {
        *l = mid;
        if (graph.subgraph_density < edge_num / sqrt(vertices[0].size() * vertices[1].size())) {
            graph.subgraph_density = edge_num / sqrt(vertices[0].size() * vertices[1].size());
            graph.vertices[0] = vertices[0];
            graph.vertices[1] = vertices[1];
        }
    } else {
        *r = mid;
    }
}