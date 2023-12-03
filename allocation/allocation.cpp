//
// Created by yy on 12/1/23.
//

#include "allocation.h"
#include <algorithm>
#include <cmath>

void Allocation::flowExactAllocation(Graph &graph, FlowNetwork &flow, double ratio, double l, double r) {
    ui n = graph.getVerticesCount();
    ui m = graph.getEdgesCount();
    auto in_degrees = graph.getInDegrees();
    auto out_degrees = graph.getOutDegrees();
    double ratio_sqrt = sqrt(ratio);
    double mid = (l + r) / 2;

    flow = FlowNetwork(2 * n + 2);
    VertexID s = 0, t = 2 * n + 1;
    for (int i = 1; i <= n; i++) {
        if (out_degrees[i - 1] > 0) {
            flow.addEdge(s, i, m);
            flow.addEdge(i, t, m + mid / ratio_sqrt);
        }
        if (in_degrees[i - 1] > 0) {
            flow.addEdge(s, i + n, m);
            flow.addEdge(i + n, t, m + ratio_sqrt * mid - 2 * in_degrees[i - 1]);
        }
        for (auto &v: graph.getOutNeighbors(i - 1)) {
            flow.addEdge(v + 1 + n, i, 2);
        }
    }
    flow.getMaxFlow(s, t);
}

//void Allocation::flowExactAllocation(Graph &graph, FlowNetwork &flow, double ratio) {
//    double l = graph.subgraph_density;
//    double r = graph.subgraph_density_upper_bound;
//    ui n = graph.getVerticesCount();
//    ui m = graph.getEdgesCount();
//    auto in_degrees = graph.getInDegrees();
//    auto out_degrees = graph.getOutDegrees();
//    if (r == 0) {
//        ui max_out = *std::max_element(out_degrees.begin(), out_degrees.end());
//        ui max_in = *std::max_element(in_degrees.begin(), in_degrees.end());
//        r = std::max(max_in, max_out);
//        graph.subgraph_density_upper_bound = r;
//    }
//
//    double ratio_sqrt = sqrt(ratio);
//    double bias = 1 / sqrt((double) n * (n - 1)) - 1 / n;
//    double mid;
//    while (r - l > bias) {
//        mid = (l + r) / 2;
//        flow = FlowNetwork(2 * n + 2);
//        VertexID s = 0, t = 2 * n + 1;
//        for (int i = 1; i <= n; i++) {
//            if (out_degrees[i - 1] > 0) {
//                flow.addEdge(s, i, m);
//                flow.addEdge(i, t, m + mid / ratio_sqrt);
//            }
//            if (in_degrees[i - 1] > 0) {
//                flow.addEdge(s, i + n, m);
//                flow.addEdge(i + n, t, m + ratio_sqrt * mid - 2 * in_degrees[i - 1]);
//            }
//            for (auto &v: graph.getOutNeighbors(i - 1)) {
//                flow.addEdge(v + 1 + n, i, 2);
//            }
//        }
//        flow.getMaxFlow(s, t);
//        std::vector<VertexID> tmp_S, tmp_s, tmp_t;
//        flow.getMinCut(s, t, tmp_S);
//        ui edge_num = 0;
//        for (auto &v: tmp_S) {
//            if (v == 0) continue;
//            if (v <= n) {
//                tmp_s.push_back(v - 1);
//                for (auto &edge: flow.adj_[v]) {
//                    if (edge.to > n && edge.to <= 2 * n && flow.dist_[edge.to] >= n) {
//                        edge_num++;
//                    }
//                }
//            } else tmp_t.push_back(v - n - 1);
//        }
//        if (!tmp_s.empty() && !tmp_t.empty()) {
//            l = mid;
//            if (graph.subgraph_density < edge_num / sqrt(tmp_s.size() * tmp_t.size())) {
//                graph.subgraph_density = edge_num / sqrt(tmp_s.size() * tmp_t.size());
//                graph.vertices[0] = tmp_s;
//                graph.vertices[1] = tmp_t;
//            }
//        } else {
//            r = mid;
//        }
//    }
//    printf("%f\n", graph.subgraph_density);
//}


