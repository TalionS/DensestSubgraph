//
// Created by yy on 23-11-15.
//

#ifndef DENSESTSUBGRAPH_REDUCTION_H
#define DENSESTSUBGRAPH_REDUCTION_H

#include "graph.h"
#include "types.h"

#include <queue>
#include <vector>
#include <iostream>
#include <list>

class Reduction {
public:
    std::vector<ui> coreDecomposition(const Graph& graph);

};

//std::vector<ui> Reduction<IsDirected>::coreDecomposition(const Graph<IsDirected>& graph){
//    const ui n = graph.getVerticesCount();
//    std::vector<ui> degree = graph.getDegrees();
//    std::vector<ui> core(n, 0);
//    std::vector<VertexID> position(n);
//
//    ui max_degree = *std::max_element(degree.begin(), degree.end());
//    std::vector<std::list<VertexID>> bins(max_degree + 1);
//
//    for (VertexID i = 0; i < n; ++i) {
//        bins[degree[i]].push_back(i);
//    }
//
//    ui idx = 0;
//    for (auto& bin : bins) {
//        for (VertexID v : bin) {
//            position[v] = idx;
//            core[idx++] = degree[v];
//        }
//    }
//
//    for (ui d = 0; d <= max_degree; ++d) {
//        for (auto it = bins[d].begin(); it != bins[d].end(); ++it) {
//            VertexID v = *it;
//            core[v] = d;
//            for (VertexID u : graph.getNeighbors(v)) {
//                if (degree[u] > d) {
//                    auto& current_bin = bins[degree[u]];
//                    auto it_u = std::find(current_bin.begin(), current_bin.end(), u);
//                    if (it_u != current_bin.end()) {
//                        current_bin.erase(it_u);
//                    }
//
//                    degree[u]--;
//                    bins[degree[u]].push_back(u);
//                }
//            }
//        }
//    }
//
//    return core;
//}

#endif //DENSESTSUBGRAPH_REDUCTION_H
