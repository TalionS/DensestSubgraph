//
// Created by yy on 12/1/23.
//

#include "allocation.h"
#include <algorithm>
#include <cmath>

void Allocation::flowExactAllocation(Graph &graph, Graph &x_y_core, FlowNetwork &flow, std::pair<double, double> ratio, double l, double r, bool is_dc) {
    ui n = graph.getVerticesCount();
    ui m = x_y_core.getEdgesCount();
    auto in_degrees = x_y_core.getInDegrees();
    auto out_degrees = x_y_core.getOutDegrees();
    double ratio_sqrt;
    if(is_dc)
        ratio_sqrt = sqrt((ratio.first + ratio.second) / 2);
    else
        ratio_sqrt = sqrt(ratio.first / ratio.second);
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

void Allocation::coreApproAllocation(Graph &graph, std::pair<ui, ui> &max_core_num_pair) {
    unsigned long long max_prod = 0;
    ui core_nums[2];
    xycore.xyCoreInitialization(graph);
    std::vector<ui> max_degrees(2, 0);
    for(int i = 0; i < 2; i++)
        max_degrees[i] = *std::max_element(xycore.degrees[i].begin(), xycore.degrees[i].end());
    ui cur = (max_degrees[0] < max_degrees[1]) ? 0 : 1;
    core_nums[cur] = 1;
    core_nums[1 - cur] = max_degrees[1 - cur];
    max_prod = core_nums[1 - cur];

    ui delta = xycore.getDelta(graph);
    printf("delta %d\n", delta);
    for (int i = 0; i < 2; i++){
        for (ui d = 1; d <= delta; d++){
            if (((unsigned long long) d) * max_degrees[1 - i] <= max_prod)
                continue;
            ui d_opst;
            d_opst = xycore.skyline_core_num(graph, i, d, static_cast<ui>(max_prod / d + 1));
            if (((unsigned long long) d) * d_opst > max_prod) {
                max_prod = d * d_opst;
                core_nums[i] = d;
                core_nums[1 - i] = d_opst;
            }
        }
    }

    printf("alpha %d beta %d\n", core_nums[0], core_nums[1]);
    max_core_num_pair =  std::make_pair(core_nums[0], core_nums[1]);
}

void Allocation::UndirectedflowExactAllocation(Graph &graph, FlowNetwork &flow, double l, double r) {
    ui n = graph.getVerticesCount();
    ui m = graph.getEdgesCount();
    auto degrees = graph.getDegrees();
    double mid = (l + r) / 2;
    flow = FlowNetwork(n + 2);
    VertexID s = 0, t = n + 1;
    for (int i = 1; i <= n; i++) {
        flow.addEdge(s, i, m);
        flow.addEdge(i, t, m + 2 * mid - degrees[i - 1]);
        for (auto &v: graph.getNeighbors(i - 1)) {
            flow.addEdge(i, v + 1, 1);
        }
    }
    flow.getMaxFlow(s, t);
}

void Allocation::UndirectedlpAllocation(Graph &graph, LinearProgamming &lp, ui T) {
    double learning_rate;
    for(ui t = T >> 1; t < T; t++){
        learning_rate = 2.0 / (t + 2);
        lp.Iterate(learning_rate);
    }
    /*
    for(int i = 0; i < graph.getVerticesCount(); i++){
        std::cout<<lp.r[0][i]<<" ";
    }
    puts("");
    for(int i = 0; i < graph.getEdgesCount(); i++){
        std::cout<<lp.alpha[0][i].id_first<<" "<<lp.alpha[0][i].weight_first<<" "<<lp.alpha[0][i].id_second<<" "<<lp.alpha[0][i].weight_second<<std::endl;
    }
    */
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
