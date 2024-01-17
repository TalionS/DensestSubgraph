//
// Created by yy on 12/1/23.
//

#include "allocation.h"

void Allocation::flowExactAllocation(Graph &graph, FlowNetwork &flow, std::pair<double, double> ratio, double l, double r, bool is_dc) {
    ui n = graph.getVerticesCount();
    ui m = graph.getEdgesCount();
    auto in_degrees = graph.getInDegrees();
    auto out_degrees = graph.getOutDegrees();
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

void Allocation::directedBSApproAllocation(Graph &graph, std::pair<double, double> ratio, std::vector<Heap> &heap,
                                           std::vector<std::vector<Heap::handle_type>> &handles,
                                           std::vector<std::vector<bool>> &is_peeled,
                                           ui &edges_count,
                                           std::vector<ui> &vertices_count,
                                           bool &is_init) {
    auto n = graph.getVerticesCount();
    if (!is_init) {
        is_init = true;
        edges_count = graph.getEdgesCount();
        vertices_count.resize(2, 0);
        heap.clear();
        heap.resize(2);
        handles.clear();
        handles.resize(2);
        is_peeled.clear();
        is_peeled.resize(2);
//        std::vector<std::vector<VertexID>> vert(2);
//        std::vector<std::vector<ui>> bin(2);
//        std::vector<std::vector<ui>> pos(2);
        std::vector<std::vector<ui>> degrees(2);
        degrees[0] = graph.getOutDegrees();
        degrees[1] = graph.getInDegrees();
        for (int i = 0; i < 2; i++) {
//            vert[i].clear();
//            bin[i].clear();
//            pos[i].clear();
//
//            vert[i].resize(n, 0);
//            bin[i].resize(n + 1, 0);
//            pos[i].resize(n, 0);
            handles[i].resize(n);
            is_peeled[i].resize(n, true);
//            for (VertexID v = 0; v < n; v++) {
//                ++bin[i][degrees[i][v]];
//            }
//            ui start = 0;
//            for (int d = 0; d <= n; d++) {
//                ui num = bin[i][d];
//                bin[i][d] = start;
//                start += num;
//            }
//            for (VertexID v = 0; v < n; v++) {
//                pos[i][v] = bin[i][degrees[i][v]]++;
//                vert[i][pos[i][v]] = v;
//            }
//            for (int d = n; d > 0; d--) {
//                bin[i][d] = bin[i][d - 1];
//            }
//            bin[i][0] = 0;
//            for (int j = bin[i][1]; j < n; j++) {
//                int u = vert[i][j];
//                handles[i][u] = heap[i].push(std::make_pair(degrees[i][u], u));
//                is_peeled[i][u] = false;
//            }
            for (ui u = 0; u < n; u++){
                if (degrees[i][u]){
                    handles[i][u] = heap[i].push(std::make_pair(degrees[i][u], u));
                    is_peeled[i][u] = false;
                }
            }
            for (ui u = 0; u < n; u++)
                if (!is_peeled[i][u])
                    vertices_count[i]++;
        }
        return;
    }

//    auto remove_node = [&](int cur) {
//        VertexID u = heap[cur].top().second;
//        heap[cur].pop();
//        is_peeled[cur][u] = true;
//        for (auto v : cur? graph.getInNeighbors(u): graph.getOutNeighbors(u)) {
//            if (!is_peeled[1 - cur][v]) {
//                (*handles[1 - cur][v]).first--;
//                heap[1 - cur].increase(handles[1 - cur][v]);
//                edges_count -= 1;
//            }
//        }
//    };
    ui  cur = heap[0].top().first * ratio.first <= heap[1].top().first * ratio.second? 0 : 1;
    VertexID u = heap[cur].top().second;
    heap[cur].pop();
    is_peeled[cur][u] = true;
    for (auto v : cur? graph.getInNeighbors(u): graph.getOutNeighbors(u)) {
        if (!is_peeled[1 - cur][v]) {
            (*handles[1 - cur][v]).first--;
            heap[1 - cur].increase(handles[1 - cur][v]);
            edges_count -= 1;
        }
    }
//    if (heap[0].top().first * ratio.first <= heap[1].top().first * ratio.second) {
//        remove_node(0);
//        --vertices_count[0];
//    } else {
//        remove_node(1);
//        --vertices_count[1];
//    }

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

void Allocation::directedKSApproAllocation(Graph &graph, std::vector<Heap> &heap,
                                           std::vector<std::vector<Heap::handle_type>> &handles,
                                           std::vector<std::vector<bool>> &is_peeled,
                                           ui &edges_count,
                                           bool &is_init) {
    ui n = graph.getVerticesCount();
    if (!is_init) {
        is_init = true;
        edges_count = graph.getEdgesCount();
        heap.clear();
        heap.resize(2);
        handles.clear();
        handles.resize(2);
        is_peeled.clear();
        is_peeled.resize(2);
//        std::vector<std::vector<VertexID>> vert(2);
//        std::vector<std::vector<ui>> bin(2);
//        std::vector<std::vector<ui>> pos(2);
        std::vector<std::vector<ui>> degrees(2);
        degrees[0] = graph.getOutDegrees();
        degrees[1] = graph.getInDegrees();
        for (int i = 0; i < 2; i++) {
//            vert[i].clear();
//            bin[i].clear();
//            pos[i].clear();
//
//            vert[i].resize(n, 0);
//            bin[i].resize(n + 1, 0);
//            pos[i].resize(n, 0);
            handles[i].resize(n);
            is_peeled[i].resize(n, true);
//            for (VertexID v = 0; v < n; v++) {
//                ++bin[i][degrees[i][v]];
//            }
//            ui start = 0;
//            for (int d = 0; d <= n; d++) {
//                ui num = bin[i][d];
//                bin[i][d] = start;
//                start += num;
//            }
//            for (VertexID v = 0; v < n; v++) {
//                pos[i][v] = bin[i][degrees[i][v]]++;
//                vert[i][pos[i][v]] = v;
//            }
//            for (int d = n; d > 0; d--) {
//                bin[i][d] = bin[i][d - 1];
//            }
//            bin[i][0] = 0;
//            for (int j = bin[i][1]; j < n; j++) {
//                int u = vert[i][j];
//                handles[i][u] = heap[i].push(std::make_pair(degrees[i][u], u));
//                is_peeled[i][u] = false;
//            }
            for (ui u = 0; u < n; u++){
                if (degrees[i][u]){
                    handles[i][u] = heap[i].push(std::make_pair(degrees[i][u], u));
                    is_peeled[i][u] = false;
                }
            }
        }
        return;
    }
//    auto remove_node = [&](ui cur) {
//        VertexID u = heap[cur].top().second;
//        heap[cur].pop();
//        is_peeled[cur][u] = true;
//        for (auto v : cur? graph.getInNeighbors(u): graph.getOutNeighbors(u)) {
//            if (!is_peeled[1 - cur][v]) {
//                if (1 - cur == 0 && v == 44)
//                    printf("neighbor %d\n", u);
//                (*handles[1 - cur][v]).first--;
//                heap[1 - cur].increase(handles[1 - cur][v]);
//                edges_count -= 1;
//            }
//        }
//    };
//    printf("%d\n", (*handles[1][6540]).first);
//    if (heap[0].empty())
//        return;
//    for (ui i = 0; i < 2; i++)
//        while (!heap[i].empty() && heap[i].top().first == 0)
//            remove_node(i);
//    if (heap[1].empty())
//        for (ui u = 0; u < n; u++){
//            if (!heap[0].empty())
//                printf("%d\n", (*handles[0][44]).first);
//            if (!is_peeled[0][u])
//                printf("%d, %d\n", u, (*handles[0][u]).first);
//        }
//    if (heap[0].empty() || heap[1].empty())
//        return;
    ui cur = heap[0].top().first <= heap[1].top().first? 0: 1;
    VertexID u = heap[cur].top().second;
    heap[cur].pop();
    is_peeled[cur][u] = true;
    for (auto v : cur? graph.getInNeighbors(u): graph.getOutNeighbors(u)) {
        if (!is_peeled[1 - cur][v]) {
            (*handles[1 - cur][v]).first--;
            heap[1 - cur].increase(handles[1 - cur][v]);
            edges_count -= 1;
        }
    }
//    if (heap[0].top().first <= heap[1].top().first) {
//        printf("%d\n", heap[0].top().first);
//        remove_node(0);
//    } else {
//        remove_node(1);
//    }
//    for (ui i = 0; i < 2; i++)
//        while (!heap[i].empty() && heap[i].top().first == 0)
//            remove_node(i);
}

void
Allocation::directedPMApproAllocation(Graph &graph, std::pair<double, double> ratio, double epsilon, ui &edges_count,
                                      std::vector<std::vector<VertexID>> &vertices,
                                      std::vector<std::vector<ui>> &degrees,
                                      bool &is_init) {
    ui n = graph.getVerticesCount();
//    double ratio = ratio.first / ratio.second;
    if (!is_init){
        is_init = true;
        edges_count = graph.getEdgesCount();
        degrees[0] = graph.getOutDegrees();
        degrees[1] = graph.getInDegrees();
        for (ui i = 0; i < 2; i++) {
            vertices[i].clear();
            for (ui u = 0; u < n; u++)
                if (degrees[i][u])
                    vertices[i].push_back(u);
        }
    }
//    double density = edges_count / sqrt(vertices[0].size() * vertices[1].size());
//    if (vertices[0].size() >= vertices[1].size() * ratio.first)
    std::vector<ui> cnt(2);
    for (ui i = 0; i < 2; i++)
        cnt[i] = vertices[i].size();
    ui side =  cnt[0] >= cnt[1] * ratio.first? 0: 1;
    std::vector<VertexID> tmp;
    ui edges_peeled_count = 0;
    for (auto u: vertices[side]) {
        if (degrees[side][u] > (1 + epsilon) * edges_count / cnt[side])
            tmp.push_back(u);
        else {
            for (auto v: side? graph.getInNeighbors(u): graph.getOutNeighbors(u))
                degrees[1 - side][v]--;
            edges_peeled_count += degrees[side][u];
        }
    }
    vertices[side] = tmp;
    edges_count -= edges_peeled_count;
}

void Allocation::directedCPAllocation(Graph &graph, LinearProgramming &lp, ui T, bool &is_init, std::pair<double, double> ratios, bool is_vw_appro) {
    double ratio;
    if (!is_vw_appro) {
        if (ratios.first < 1 && ratios.second > 1) {
            ratio = 1;
        } else if (ratios.second <= 1) {
            ratio = (ratios.first + ratios.second) / 2;
        } else if (ratios.first >= 1) {
            ratio = 2 / (1 / ratios.first + 1 / ratios.second);
        }
    } else{
        ratio = ratios.first / ratios.second;
    }
    if (!is_init) {
        lp.Init(graph, ratio);
        is_init = true;
    }
    double learning_rate;
//    for (ui t = T - 100; t < T; t++){
    for (ui t = T >> 1; t < T; t++){
        learning_rate = 2.0 / (t + 2);
        lp.Iterate(learning_rate, ratio);
    }
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

void Allocation::UndirectedlpAllocation(Graph &graph, LinearProgramming &lp, ui T) {
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
