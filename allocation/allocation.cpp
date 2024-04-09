//
// Created by yy on 12/1/23.
//

#include "allocation.h"

void
Allocation::flowExactAllocation(Graph &graph, FlowNetwork &flow, std::pair<double, double> ratio, double l, double r,
                                bool is_dc, bool is_map) {
    ui n = graph.getVerticesCount();
    ui m = graph.getEdgesCount();
    std::vector<ui> degrees[2];
    degrees[1] = graph.getInDegrees();
    degrees[0] = graph.getOutDegrees();
    double ratio_sqrt;
    if (is_dc) {
        double ratios;
        if (is_map) {
            if (ratio.first < 1 && ratio.second > 1) {
                ratios = 1;
            } else if (ratio.second <= 1) {
                ratios = (ratio.first + ratio.second) / 2;
            } else if (ratio.first >= 1) {
                ratios = 2 / (1 / ratio.first + 1 / ratio.second);
            }
        } else
            ratios = (ratio.first + ratio.second) / 2;
        ratio_sqrt = sqrt(ratios);
    }
    else
        ratio_sqrt = sqrt(ratio.first / ratio.second);
    double mid = (l + r) / 2;
    flow = FlowNetwork(2 * n + 2);
    VertexID s = 0, t = 2 * n + 1;
    for (int i = 1; i <= n; i++) {
        if (degrees[0][i - 1] > 0) {
            flow.addEdge(s, i, m);
            flow.addEdge(i, t, m + mid / ratio_sqrt);
        }
        if (degrees[1][i - 1] > 0) {
            flow.addEdge(s, i + n, m);
            flow.addEdge(i + n, t, m + ratio_sqrt * mid - 2 * degrees[1][i - 1]);
        }
        for (auto &v: graph.getOutNeighbors(i - 1)) {
            flow.addEdge(v + 1 + n, i, 2);
        }
    }
    flow.getMaxFlow(s, t);
//    std::vector<ui> map[2];
//    ui tmp = 1;
//    for (ui i = 0; i < 2; i++){
//        map[i].resize(n);
//        for (VertexID u = 0; u < n; u++){
//            if (degrees[i][u])
//                map[i][u] = tmp++;
//        }
//    }
//    VertexID s = 0, t = tmp;
//    flow = FlowNetwork(t + 1);
//    for (int i = 0; i < n; i++) {
//        if (degrees[0][i] > 0) {
//            flow.addEdge(s, map[0][i], m);
//            flow.addEdge(map[0][i], t, m + mid / ratio_sqrt);
//        }
//        if (degrees[1][i] > 0) {
//            flow.addEdge(s, map[1][i], m);
//            flow.addEdge(map[1][i], t, m + ratio_sqrt * mid - 2 * degrees[1][i]);
//        }
//        for (auto &v: graph.getOutNeighbors(i)) {
//            flow.addEdge(map[1][v], map[0][i], 2);
//        }
//    }
//    flow.getMaxFlow(s, t);
}

void Allocation::directedBSApproAllocation(Graph &graph, std::pair<double, double> ratio, std::vector<Heap> &heap,
                                           std::vector<std::vector<Heap::handle_type>> &handles,
                                           std::vector<std::vector<bool>> &is_peeled,
                                           ui &edges_count,
                                           bool &is_init) {
    auto n = graph.getVerticesCount();
    if (!is_init) {
        is_init = true;
        edges_count = graph.getEdgesCount();
//        vertices_count.resize(2, 0);
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
            handles[i].resize(n);
            is_peeled[i].resize(n, true);
            for (ui u = 0; u < n; u++) {
                if (degrees[i][u]) {
                    handles[i][u] = heap[i].push(std::make_pair(-degrees[i][u], u));
                    is_peeled[i][u] = false;
                }
            }
        }
        return;
    }

    ui cur = heap[0].top().first * ratio.first >= heap[1].top().first * ratio.second ? 0 : 1;
    VertexID u = heap[cur].top().second;
    heap[cur].pop();
    is_peeled[cur][u] = true;
    for (auto v: cur ? graph.getInNeighbors(u) : graph.getOutNeighbors(u)) {
        if (!is_peeled[1 - cur][v]) {
            (*handles[1 - cur][v]).first++;
            heap[1 - cur].increase(handles[1 - cur][v]);
            edges_count -= 1;
        }
    }

}

void Allocation::wCoreApproAllocation(Graph &graph, WCore &w_core, std::pair<ui, ui> &max_core_num_pair) {
    w_core.getMaxCNPair(graph, max_core_num_pair);
}

void Allocation::xyCoreApproAllocation(Graph &graph, std::pair<ui, ui> &max_core_num_pair) {
    unsigned long long max_prod = 0;
    ui core_nums[2];
    xycore.xyCoreInitialization(graph, true);
    std::vector<ui> max_degrees(2, 0);
    for (int i = 0; i < 2; i++)
        max_degrees[i] = *std::max_element(xycore.degrees[i].begin(), xycore.degrees[i].end());
    ui cur = (max_degrees[0] < max_degrees[1]) ? 0 : 1;
    core_nums[cur] = 1;
    core_nums[1 - cur] = max_degrees[1 - cur];
    max_prod = core_nums[1 - cur];
//    clock_t begin = clock();
    ui delta = xycore.getDelta(graph);
//    clock_t end = clock();
//    printf("time: %f\n", (double) (end - begin) / CLOCKS_PER_SEC);
//    printf("delta %d\n", delta);
//    begin = clock();
    for (int i = 0; i < 2; i++) {
        for (ui d = 1; d <= delta; d++) {
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
//    end = clock();
//    printf("time: %f\n", (double) (end - begin) / CLOCKS_PER_SEC);
//    printf("alpha %d beta %d\n", core_nums[0], core_nums[1]);
    max_core_num_pair = std::make_pair(core_nums[0], core_nums[1]);
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
        std::vector<std::vector<ui>> degrees(2);
        degrees[0] = graph.getOutDegrees();
        degrees[1] = graph.getInDegrees();
        for (int i = 0; i < 2; i++) {
            handles[i].resize(n);
            is_peeled[i].resize(n, true);
            for (ui u = 0; u < n; u++) {
                if (degrees[i][u]) {
                    handles[i][u] = heap[i].push(std::make_pair(-degrees[i][u], u));
                    is_peeled[i][u] = false;
                }
            }
        }
        return;
    }
    ui cur = heap[0].top().first >= heap[1].top().first ? 0 : 1;
    VertexID u = heap[cur].top().second;
    heap[cur].pop();
    is_peeled[cur][u] = true;
    for (auto v: cur ? graph.getInNeighbors(u) : graph.getOutNeighbors(u)) {
        if (!is_peeled[1 - cur][v]) {
            (*handles[1 - cur][v]).first++;
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

void Allocation::directedFixedKSApproAllocation(Graph &graph, std::pair<double, double> ratio, ui cur,
                                                std::vector<Heap> &heap,
                                                std::vector<std::vector<Heap::handle_type>> &handles,
                                                std::vector<std::vector<bool>> &is_peeled, ui &edges_count) {
    ui n = graph.getVerticesCount();
    edges_count = graph.getEdgesCount();
    heap.clear();
    heap.resize(2);
    handles.clear();
    handles.resize(2);
    is_peeled.clear();
    is_peeled.resize(2);
    std::vector<std::vector<ui>> degrees(2);
    degrees[0] = graph.getOutDegrees();
    degrees[1] = graph.getInDegrees();
    for (int i = 0; i < 2; i++) {
        handles[i].resize(n);
        is_peeled[i].resize(n, true);
        for (ui u = 0; u < n; u++) {
            handles[i][u] = heap[i].push(std::make_pair(-degrees[i][u], u));
            is_peeled[i][u] = false;
        }
    }

    std::vector<ui> cnt(2, n);
    std::vector<ui> required(2);
    required[0] = (ui) ratio.first;
    required[1] = (ui) ratio.second;

    while (required[cur] < cnt[cur]) {
//        printf("vertex: %d, degrees: %d\n", heap[cur].top().second, heap[cur].top().first);
//        heap[cur].pop();
//        continue;
        VertexID u = heap[cur].top().second;
        heap[cur].pop();
        is_peeled[cur][u] = true;
        for (auto v: cur ? graph.getInNeighbors(u) : graph.getOutNeighbors(u)) {
            if (!is_peeled[1 - cur][v]) {
                (*handles[1 - cur][v]).first++;
                heap[1 - cur].increase(handles[1 - cur][v]);
                edges_count -= 1;
            }
        }
        cnt[cur]--;
    }

    while (required[1 - cur] < cnt[1 - cur]) {
        VertexID u = heap[1 - cur].top().second;
        heap[1 - cur].pop();
        is_peeled[1 - cur][u] = true;
        for (auto v: (1 - cur) ? graph.getInNeighbors(u) : graph.getOutNeighbors(u)) {
            if (!is_peeled[cur][v]) {
                (*handles[cur][v]).first++;
                heap[cur].increase(handles[cur][v]);
                edges_count -= 1;
            }
        }
        cnt[1 - cur]--;
    }
}

void
Allocation::directedPMApproAllocation(Graph &graph, std::pair<double, double> ratio, double epsilon, ui &edges_count,
                                      std::vector<std::vector<VertexID>> &vertices,
                                      std::vector<std::vector<ui>> &degrees,
                                      bool &is_init) {
    ui n = graph.getVerticesCount();
//    double ratio = ratio.first / ratio.second;
    if (!is_init) {
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
    ui side = cnt[0] >= cnt[1] * ratio.first ? 0 : 1;
    std::vector<VertexID> tmp;
    ui edges_peeled_count = 0;
    for (auto u: vertices[side]) {
        if (degrees[side][u] > (1 + epsilon) * edges_count / cnt[side])
            tmp.push_back(u);
        else {
            for (auto v: side ? graph.getInNeighbors(u) : graph.getOutNeighbors(u))
                degrees[1 - side][v]--;
            edges_peeled_count += degrees[side][u];
        }
    }
    vertices[side] = tmp;
    edges_count -= edges_peeled_count;
}

void
Allocation::directedCPAllocation(Graph &graph, LinearProgramming &lp, ui &iter_num, bool &is_init,
                                 std::pair<double, double> ratios, bool is_synchronous, bool is_exp, bool is_map) {
    double ratio;
    if (is_map) {
        if (ratios.first < 1 && ratios.second > 1) {
            ratio = 1;
        } else if (ratios.second <= 1) {
            ratio = (ratios.first + ratios.second) / 2;
        } else if (ratios.first >= 1) {
            ratio = 2 / (1 / ratios.first + 1 / ratios.second);
        }
    } else
        ratio = (ratios.first + ratios.second) / 2;
    if (!is_init) {
        lp.Init(graph, ratio);
        is_init = true;
    }
    double learning_rate;
//    for (ui t = T - 100; t < T; t++){
    ui cur_iter_num = lp.cur_iter_num;
    if (is_exp)
        iter_num = cur_iter_num? cur_iter_num: 1;
    for (ui t = cur_iter_num; t < cur_iter_num + iter_num; t++) {
        learning_rate = 2.0 / (t + 2);
        lp.Iterate(learning_rate, ratio, is_synchronous);
    }
//    printf("%d\n", lp.cur_iter_num);
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

void Allocation::UndirectedlpAllocation(Graph &graph, LinearProgramming &lp, ui T, bool is_synchronous) {
    double learning_rate;
    for (ui t = T >> 1; t < T; t++) {
        learning_rate = 2.0 / (t + 2);
        lp.Iterate(learning_rate, 0, is_synchronous);
    }
    /*
    for(int i = 0; i < graph.getVerticesCount(); i++){
        std::cout<<lp.r[0][i]<<" ";
    }
    puts("");
    for(int i = 0; i < graph.getEdgesCount(); i++){
        std::cout<<lp.alpha[i].id_first<<" "<<lp.alpha[i].weight_first<<" "<<lp.alpha[i].id_second<<" "<<lp.alpha[i].weight_second<<std::endl;
    }
    */
}

void Allocation::UndirectedFistaAllocation(Graph &graph, LinearProgramming &lp, ui T, bool is_synchronous) {
    double max_deg = graph.getMaxdeg();
    for (ui t = T >> 1; t < T; t++) {
        lp.FistaIterate(1.0 / 2 / max_deg, t + 1, 0, is_synchronous);
    }
}

void Allocation::UndirectedMWUAllocation(Graph &graph, LinearProgramming &lp, ui T, bool is_synchronous) {
    for (ui t = T >> 1; t < T; t++) {
        lp.MWUIterate(t + 1, is_synchronous);
    }
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
void Allocation::UndirectedGreedyAllocation(Graph &graph){
    std::vector<ui> deg = graph.getDegrees();
    ui num_vertex = graph.getVerticesCount();
    ui num_edge = graph.getEdgesCount();
    ui n = graph.getVerticesCount();
    double opt = 1.0 * num_edge / num_vertex;
    std::vector<ui> *bucket = new std::vector<ui>[n]();
    std::vector<ui> del;
    del.resize(n);
    for(ui i = 0; i < n; i++){
        bucket[deg[i]].push_back(i);
        del[i] = 0;
    }
    for(ui i = 0; i < n; i++){
        for(int j = 0; j < bucket[i].size(); j++){
            ui u = bucket[i][j];
            if(del[u] == 1) continue;
            del[u] = 1;
            num_vertex--;
            for(auto v : graph.getNeighbors(u)){
                if(del[v] == 1) continue;
                deg[v]--;
                num_edge--;
                bucket[deg[v]].push_back(v);
            }
        }
        opt = std::max(opt, 1.0 * num_edge / num_vertex);
    }
    std::cout<<opt<<std::endl;
}
void Allocation::UndirectedCoreAppAllocation(Graph &graph, CoreApp &ca) {
    ui kmax = 0, kmin = ca.nodes_count + 1;
    ui num_edge = 0;
    while(ca.pos < ca.size && ca.pos != ca.nodes_count){
        ui maxx = 0;
        for(ui i = 0; i < ca.nodes_count; i++){
            if(ca.selected[i]) continue;
            maxx = std::max(maxx, ca.deg[i]);
        }
        for(ui i = 0; i < ca.nodes_count; i++){
            if(ca.deg[i] >= maxx && ca.selected[i] == 0){
                ca.selected[i] = 1;
                ca.id[ca.pos++] = i;
                for(auto v: graph.getNeighbors(i)){
                    if(ca.selected[v]){
                        ca.new_edge[v].push_back(i);
                        ca.new_edge[i].push_back(v);
                    }
                }
            }
        }
        std::cout<<maxx<<"pos = "<<ca.pos<<std::endl;
    }
    for(ui i = 0; i < ca.pos; i++){
        ui u = ca.id[i];
        ca.w[u] = 0;
        for(auto v : ca.new_edge[u]){
            if(ca.selected[v] == 0) continue;
            ca.w[u]++;
            if(v < u) num_edge++;
        }
        kmax = std::max(kmax, ca.w[u]);
        kmin = std::min(kmin, ca.w[u]);
    }
    ui num_vertex = ca.pos;
    std::vector<ui> *bucket = new std::vector<ui>[kmax + 1];
    std::vector<ui> del;
    del.resize(ca.nodes_count);
    for(ui i = 0; i < ca.pos; i++){
        ui u = ca.id[i];
        bucket[ca.w[u]].push_back(u);
        del[u] = 0;
    }
    ui pos = 0;
    ui k = std::max(kmin, ca.k + 1);
    while(k <= kmax && num_vertex){ 
        for(ui i = pos; i < k; i++){
            for(int j = 0; j < bucket[i].size(); j++){
                ui u = bucket[i][j];
                if(del[u] == 1) continue;
                del[u] = 1;
                num_vertex--;
                for(auto v : ca.new_edge[u]){
                    if(del[v] == 1) continue;
                    ca.w[v]--;
                    num_edge--;
                    bucket[ca.w[v]].push_back(v);
                }
            }
        }
        pos = k;
        if(num_vertex){
            if(k > ca.k){
                ca.k = k;
                if(1.0 * num_edge / num_vertex > ca.opt) ca.opt = 1.0 * num_edge / num_vertex;
            }
            k = k + 1;
        }
    }
    std::cout<<"pos = "<<ca.pos<<std::endl;
}

void Allocation::UndirectedGreedyppAllocation(Graph &graph, Greedy &gr){
    std::vector<ui> deg = graph.getDegrees();
    ui num_vertex = graph.getVerticesCount();
    ui num_edge = graph.getEdgesCount();
    ui n = graph.getVerticesCount();
    double opt = 1.0 * num_edge / num_vertex;
    std::vector<ui> del;
    del.resize(n);
    BinaryHeap Heap;
    Heap.resize(n);
    for(ui i = 0; i < n; i++){
        gr.h[i] = gr.h[i] + deg[i];
        del[i] = 0;
    }
    Heap.init(gr.h);
    for(ui i = 0; i < n; i++){
        HeapElement tmp = Heap.pop();
        ui u = tmp.id;
        gr.h[u] = tmp.val;
        del[u] = 1;
        num_vertex--;
        for(auto v : graph.getNeighbors(u)){
            if(del[v]) continue;
            Heap.modify(v);
            num_edge--;
        }
        if(num_edge){
            opt = std::max(opt, 1.0 * num_edge / num_vertex);
        }
    }
    graph.subgraph_density_lower_bound = std::max(graph.subgraph_density_lower_bound, opt);
}

void Allocation::UndirectedFlowAppAllocation(Graph &graph){

    double mid = (graph.subgraph_density_lower_bound + graph.subgraph_density_upper_bound) / 2;
    double eps = (mid - graph.subgraph_density_lower_bound) / (2 * graph.subgraph_density_lower_bound);
    ui n = graph.getVerticesCount();
    ui m = graph.getEdgesCount();
    ui T = n + m + 1;
    Dinic dinic;
    dinic.Init(T + 1);
    ui pos = n;
    for(ui i = 0; i < n; i++){
        dinic.addEdge(i + 1, T, mid);
        for(auto v : graph.getNeighbors(i)){
            if(i > v) continue;
            pos++;
            dinic.addEdge(0, pos, 1);
            dinic.addEdge(pos, i + 1, 1);
            dinic.addEdge(pos, v + 1, 1);
        }
    }
    double iter = ceil(2*log(m)/eps) + 2;
    double total = dinic.solve(iter);
    if(total + 0.000005 < m){
        /*
        std::vector<bool> used;
        used.resize(n);
        double epss = 1e-6;
        for(ui i = 0; i < dinic.e[T].size(); i++){
            ui v = dinic.e[T][i].to;
            if(dinic.e[T][i].flow < mid) used[v - 1] = 0;
            used[v - 1] = 1;
        }
        double num_vertex = 0;
        double num_edge = 0;
        for(ui i = 0; i < n ; i++){
            if(used[i]) num_vertex++;
            else continue;
            for(auto v : graph.getNeighbors(i)){
                if(used[v]) num_edge++;
            }
        }
        */
        graph.subgraph_density_lower_bound =  (graph.subgraph_density_lower_bound + mid) / 2;
    }
    else{
        graph.subgraph_density_upper_bound = mid;
    }
}