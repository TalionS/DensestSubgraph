//
// Created by yy on 11/23/23.
//
#include "reduction.h"
#include "lp.h"
#include "heap.h"
void Reduction::coreDecomposition(const Graph &graph) {
    ui vertices_count = graph.getVerticesCount();
    std::vector<ui> degrees = graph.getDegrees();
    core.resize(vertices_count, 0);
    std::vector<ui> rid(vertices_count, 0);
    std::vector<ui> id(vertices_count, 0);
    for (ui i = 0; i < vertices_count; i++) ++id[degrees[i]];
    for (ui i = 1; i < vertices_count; i++) id[i] += id[i - 1];

    for (ui i = 0; i < vertices_count; i++) rid[i] = --id[degrees[i]];
    for (ui i = 0; i < vertices_count; i++) id[rid[i]] = i;

    std::vector<ui> degree_start(vertices_count + 1);
    for (ui i = 0, j = 0; i <= vertices_count; i++) {
        while (j < vertices_count && degrees[id[j]] < i) ++j;
        degree_start[i] = j;
    }

    ui max_core = 0;
    for (ui i = 0; i < vertices_count; i++) {
        ui u = id[i];
//        assert(degree_start[deg[u]] == i);
        if (degrees[u] > max_core) max_core = degrees[u];
        core[u] = max_core;

        ++degree_start[degrees[u]];
        if (degrees[u] == 0) continue;

        degree_start[degrees[u] - 1] = degree_start[degrees[u]];
        for (ui j: graph.getNeighbors(u)) {
            if (rid[j] > i) {
                ui pos1 = degree_start[degrees[j]], pos2 = rid[j];
                std::swap(id[pos1], id[pos2]);
                rid[id[pos1]] = pos1;
                rid[id[pos2]] = pos2;
                ++degree_start[degrees[j]];
                --degrees[j];
            }
        }
    }
}

void Reduction::xyCoreReduction(Graph &graph, Graph &x_y_core, std::pair<double, double> ratios, double &l, double &r,
                                bool &is_init, bool is_dc, bool is_divide_by_number, bool is_exact, bool is_map,
                                bool is_res, ui res_width, bool is_copy) {
    if (!is_init) {
        is_init = true;
        xycore.xyCoreInitialization(graph, true);
    }
    ui x, y;
    if (!is_dc) {
        double ratio_sqrt = sqrt(ratios.first / ratios.second);
        x = std::max(static_cast<int>(ceil(l / 2 / ratio_sqrt)), 1);
        y = std::max(static_cast<int>(ceil(ratio_sqrt * l / 2)), 1);
    } else {
        double ratio;
        if (is_divide_by_number) {
            if (ratios.first < 1 && ratios.second > 1) {
                ratio = 1;
            } else if (ratios.second <= 1) {
                ratio = (ratios.first + ratios.second) / 2;
            } else if (ratios.first >= 1) {
                ratio = 2 / (1 / ratios.first + 1 / ratios.second);
            }
        } else
            ratio = (ratios.first + ratios.second) / 2;
//        ratios.first = std::max(ratios.first, ratio / 2);
//        double ratio_right_sqrt = sqrt(std::min(ratios.second, ratio * 2));
//        double ratio_left_sqrt = sqrt(std::max(ratios.first, ratio / 2));
        double ratio_right_sqrt, ratio_left_sqrt;
        if (!is_res){
            ratio_right_sqrt = sqrt(ratios.second);
            ratio_left_sqrt = sqrt(ratios.first);
        } else {
            ratio_right_sqrt = sqrt(std::min(ratios.second, ratio * res_width));
            ratio_left_sqrt = sqrt(std::max(ratios.first, ratio / res_width));
        }
        x = std::max(static_cast<int>(ceil(l / 2 / ratio_right_sqrt)), 1);
        y = std::max(static_cast<int>(ceil(ratio_left_sqrt * l / 2)), 1);
//        printf("x: %d, y: %d\n", x, y);
//        x = std::min(x, xycore.max_degrees[0]);
//        y = std::min(y, xycore.max_degrees[1]);
    }

    xycore.generateXYCore(graph, x_y_core, x, y, is_exact, is_map, is_copy);
}
void Reduction::stableSetReduction(Graph &graph, LinearProgramming &lp,
                                   std::vector<std::pair<VertexID, VertexID>> &edges, bool stable_set_reduction, bool is_map) {
    if (!stable_set_reduction)
        return;
//    Graph stable_subgraph = Graph(true, graph.getVerticesCount());
    if (is_map){
        auto n = graph.getVerticesCount();
        std::vector<bool> selected[2];
        for (ui i = 0; i < 2; i++)
            selected[i].resize(n, false);
        for (auto &edge: edges) {
//        stable_subgraph.addDirectedEdge(edge.first, edge.second);
            selected[0][edge.first] = true;
            selected[1][edge.second] = true;
        }
        ui vertex_id = 0;
        std::vector<VertexID> map(n);
        std::vector<VertexID> reverse_map;
//        if (graph.is_mapped) {
//            for (ui u = 0; u < n; u++) {
//                if (selected[0][u] || selected[1][u]) {
//                    reverse_map.emplace_back(graph.map[u]);
//                    map[u] = vertex_id++;
//                }
//            }
//        } else {
            for (ui u = 0; u < n; u++) {
                if (selected[0][u] || selected[1][u]) {
                    reverse_map.emplace_back(u);
                    map[u] = vertex_id++;
                }
            }
//        }
        Graph stable_subgraph = Graph(true, vertex_id);
        stable_subgraph.is_mapped = true;
//        stable_subgraph.map = reverse_map;
        for (auto &edge: edges)
            stable_subgraph.addDirectedEdge(map[edge.first], map[edge.second]);
        ui i = 0, j = lp.edges_count_ - 1;
        while (i != j) {
            while (i < j && selected[0][lp.alpha[i].id_first] && selected[1][lp.alpha[i].id_second])
//            lp.alpha[i].id_first = map[lp.alpha[i].id_first];
//            lp.alpha[i].id_second = map[lp.alpha[i].id_second];
                i++;
            while (j > i && (!selected[0][lp.alpha[j].id_first] || !selected[1][lp.alpha[j].id_second]))
                j--;
//        printf("%d, %d\n", i, j);
            std::swap(lp.alpha[i], lp.alpha[j]);

        }
        std::vector<std::vector<double>> r(2);
        for (i = 0; i < 2; i++) {
            for (ui u = 0; u < vertex_id; u++)
                r[i].emplace_back(lp.r[i][reverse_map[u]]);
        }
        if (graph.is_mapped) {
            for (auto &u: reverse_map)
                u = graph.map[u];
        }
        stable_subgraph.map = reverse_map;

        lp.edges_count_ = edges.size();
        lp.nodes_count_ = vertex_id;
        lp.alpha.resize(lp.edges_count_);
        lp.r = r;
        for (i = 0; i < lp.edges_count_; i++) {
            lp.alpha[i].id_first = map[lp.alpha[i].id_first];
            lp.alpha[i].id_second = map[lp.alpha[i].id_second];
        }
//    lp.sort(stable_subgraph);
//    for (ui i = 0; i < lp.edges_count_; i++) {
//        if (!lp.alpha[i].is_selected)
//            continue;
//        if (!selected[0][lp.alpha[i].id_first] || !selected[1][lp.alpha[i].id_second]) {
//            lp.alpha[i].is_selected = false;
//        }
//    }
        graph = stable_subgraph;
//        printf("%d\n", graph.getEdgesCount());
    } else {
        auto n = graph.getVerticesCount();
        std::vector<bool> selected[2];
        for (ui i = 0; i < 2; i++)
            selected[i].resize(n, false);
        for (auto &edge: edges) {
//        stable_subgraph.addDirectedEdge(edge.first, edge.second);
            selected[0][edge.first] = true;
            selected[1][edge.second] = true;
        }
        Graph stable_subgraph = Graph(true, n);
        for (auto &edge: edges)
            stable_subgraph.addDirectedEdge(edge.first, edge.second);
//    std::vector<std::vector<ui>> deg(2);
//    deg[0] = stable_subgraph.getOutDegrees();
//    deg[1] = stable_subgraph.getInDegrees();
//    for (VertexID u = 0; u < n; u++) {
//        if (deg[0][u])
//            stable_subgraph.vertex_ids[0].emplace_back(u);
//        if (deg[1][u])
//            stable_subgraph.vertex_ids[1].emplace_back(u);
//    }
        ui i = 0, j = lp.edges_count_ - 1;
        while (i != j) {
            while (i < j && selected[0][lp.alpha[i].id_first] && selected[1][lp.alpha[i].id_second])
//            lp.alpha[i].id_first = map[lp.alpha[i].id_first];
//            lp.alpha[i].id_second = map[lp.alpha[i].id_second];
                i++;
            while (j > i && (!selected[0][lp.alpha[j].id_first] || !selected[1][lp.alpha[j].id_second]))
                j--;
//        printf("%d, %d\n", i, j);
            std::swap(lp.alpha[i], lp.alpha[j]);

        }
        lp.edges_count_ = edges.size();
//    lp.nodes_count_ = vertex_id;
        lp.alpha.resize(lp.edges_count_);
//    lp.r = r;

//    lp.sort(stable_subgraph);
//    for (ui i = 0; i < lp.edges_count_; i++) {
//        if (!lp.alpha[i].is_selected)
//            continue;
//        if (!selected[0][lp.alpha[i].id_first] || !selected[1][lp.alpha[i].id_second]) {
//            lp.alpha[i].is_selected = false;
//        }
//    }
        graph = stable_subgraph;
    }
}
//void Reduction::stableSetReduction(Graph &graph, LinearProgramming &lp,
//                                   std::vector<std::pair<VertexID, VertexID>> &edges, bool stable_set_reduction) {
//    if (!stable_set_reduction)
//        return;
////    Graph stable_subgraph = Graph(true, graph.getVerticesCount());
//    auto n = graph.getVerticesCount();
//    std::vector<bool> selected[2];
//    for (ui i = 0; i < 2; i++)
//        selected[i].resize(n, false);
//    for (auto &edge: edges) {
////        stable_subgraph.addDirectedEdge(edge.first, edge.second);
//        selected[0][edge.first] = true;
//        selected[1][edge.second] = true;
//    }
//    ui vertex_id = 0;
//    std::vector<VertexID> map(n);
//    std::vector<VertexID> reverse_map(n);
//    for (ui u = 0; u < n; u++) {
//        if (selected[0][u] || selected[1][u]) {
//            reverse_map[vertex_id] = graph.map[u];
//            map[u] = vertex_id++;
//        }
//    }
//    reverse_map.resize(vertex_id);
//    Graph stable_subgraph = Graph(true, vertex_id);
//    stable_subgraph.map = reverse_map;
//    for (auto &edge: edges)
//        stable_subgraph.addDirectedEdge(map[edge.first], map[edge.second]);
//    ui i = 0, j = lp.edges_count_ - 1;
//    while (i != j) {
//        while (i < j && selected[0][lp.alpha[i].id_first] && selected[1][lp.alpha[i].id_second])
////            lp.alpha[i].id_first = map[lp.alpha[i].id_first];
////            lp.alpha[i].id_second = map[lp.alpha[i].id_second];
//            i++;
//        while (j > i && (!selected[0][lp.alpha[j].id_first] || !selected[1][lp.alpha[j].id_second]))
//            j--;
////        printf("%d, %d\n", i, j);
//        std::swap(lp.alpha[i], lp.alpha[j]);
//
//    }
//    std::vector<std::vector<double>> r(2);
//    for (i = 0; i < 2; i++) {
//        for (ui u = 0; u < vertex_id; u++)
//            r[i].emplace_back(lp.r[i][reverse_map[u]]);
//    }
//    lp.edges_count_ = edges.size();
//    lp.nodes_count_ = vertex_id;
//    lp.alpha.resize(lp.edges_count_);
//    lp.r = r;
//    for (ui i = 0; i < lp.edges_count_; i++) {
//        lp.alpha[i].id_first = map[lp.alpha[i].id_first];
//        lp.alpha[i].id_second = map[lp.alpha[i].id_second];
//    }
////    lp.sort(stable_subgraph);
////    for (ui i = 0; i < lp.edges_count_; i++) {
////        if (!lp.alpha[i].is_selected)
////            continue;
////        if (!selected[0][lp.alpha[i].id_first] || !selected[1][lp.alpha[i].id_second]) {
////            lp.alpha[i].is_selected = false;
////        }
////    }
//    graph = stable_subgraph;
//}

void Reduction::wCoreReduction(Graph &graph, Graph &subgraph, WCore &w_core) {
    w_core.generateMaxWCore(graph, subgraph);
}

void Reduction::kCoreReduction(Graph &graph, double &l, double &r) {
    ui vertices_count = graph.getVerticesCount();
    std::vector<ui> degrees = graph.getDegrees();
    std::vector<ui> newid(vertices_count, 0);
    std::queue<VertexID> q;
    for (int i = 0; i < vertices_count; i++) {
        if (degrees[i] < l) q.push(i), newid[i] = vertices_count;
    }
    while (q.size()) {
        int u = q.front();
        q.pop();
        for (ui v: graph.getNeighbors(u)) {
            if (newid[v] == vertices_count) continue;
            if ((--degrees[v]) < l) {
                q.push(v);
                newid[v] = vertices_count;
            }
        }
    }
    ui new_vertices_count = 0;
    for (int i = 0; i < vertices_count; i++) {
        if (newid[i] != vertices_count) newid[i] = new_vertices_count++;
    }
    Graph kcore = Graph(false, new_vertices_count);
    for (int i = 0; i < vertices_count; i++) {
        if (newid[i] == vertices_count) continue;
        for (ui j: graph.getNeighbors(i)) {
            if (newid[j] == vertices_count) continue;
            if (newid[i] < newid[j]) kcore.addUndirectedEdge(newid[i], newid[j]);
        }
    }
    kcore.subgraph_density_lower_bound = l;
    kcore.subgraph_density_upper_bound = r;
    graph = kcore;
    auto dgrees = graph.getDegrees();
}

void Reduction::UndirectedkCoreReduction(Graph &graph, LinearProgramming &lp, bool weighted){
    if(!weighted){
        ui vertices_count_ = graph.getVerticesCount();
        ui bound = (ui) (graph.subgraph_density_lower_bound);
        std::vector<ui> *bucket = new std::vector<ui>[graph.getVerticesCount()];
        for(ui i = 0 ; i < vertices_count_; i++) bucket[i].clear();
        std::vector<ui> deg = graph.getDegrees();
        std::vector<ui> id;
        id.resize(vertices_count_);
        for(ui i = 0; i < vertices_count_; i++){
            bucket[deg[i]].push_back(i);
            id[i] = vertices_count_;
        }
        for(ui i = 0; i < vertices_count_; i++){
            for(ui j = 0; j < bucket[i].size(); j++){
                ui u = bucket[i][j];
                if(id[u] != vertices_count_) continue;
                id[u] = i;
                for(auto v : graph.getNeighbors(u)){
                    if(id[v] != vertices_count_) continue;
                    bucket[--deg[v]].push_back(v);
                }
            }
        }
        ui new_vertices_count_ = 0;
        for(ui i = 0; i < vertices_count_; i++){
            if(id[i] >= bound){
                id[i] = bound + new_vertices_count_;
                new_vertices_count_++;
            }        
        }
        Graph copy(0, new_vertices_count_);
        copy.subgraph_density = graph.subgraph_density;
        copy.weight_.resize(new_vertices_count_);
        for(ui i = 0; i < new_vertices_count_; i++) copy.weight_[i] = 1;
        for(ui i = 0; i < vertices_count_; i++){
            if(id[i] < bound) continue;
            for(auto v: graph.getNeighbors(i)){
                if(id[v] < bound || id[i] < id[v]) continue;
                copy.addUndirectedEdge(id[i] - bound, id[v] - bound);
            }
        }
        std::vector<Alpha> new_alpha;
        for(ui i = 0; i < graph.getEdgesCount(); i++){
            if(id[lp.alpha[i].id_first] >= bound && id[lp.alpha[i].id_second] >= bound){
                lp.alpha[i].id_first = id[lp.alpha[i].id_first] - bound;
                lp.alpha[i].id_second= id[lp.alpha[i].id_second] - bound;
                new_alpha.push_back(lp.alpha[i]);
            }
        }
        lp.alpha = new_alpha;
        if(lp.type_){
            std::vector<Alpha> new_beta;
            for(ui i = 0; i < graph.getEdgesCount(); i++){
                if(id[lp.beta[i].id_first] >= bound && id[lp.beta[i].id_second] >= bound){
                    lp.beta[i].id_first = id[lp.beta[i].id_first] - bound;
                    lp.beta[i].id_second= id[lp.beta[i].id_second] - bound;
                    new_beta.push_back(lp.beta[i]);
                }
            }
            lp.beta = new_beta;
        }
        lp.nodes_count_ = copy.getVerticesCount();
        lp.edges_count_ = copy.getEdgesCount();
        for(ui i = 0; i < lp.nodes_count_; i++) lp.r[0][i] = 0;
        for(ui i = 0; i < lp.edges_count_; i++){
            lp.r[0][lp.alpha[i].id_first] += lp.alpha[i].weight_first;
            lp.r[0][lp.alpha[i].id_second] += lp.alpha[i].weight_second;
        }
        lp.weight = copy.weight_;
        graph = copy;
    }
    else{
        BinaryHeap heap;
        ui vertices_count_ = graph.getVerticesCount(); 
        heap.resize(vertices_count_);
        std::vector<ui> deg = graph.getDegrees();
        std::vector<double> h;
        h.resize(vertices_count_);
        std::vector<ui> id;
        id.resize(vertices_count_);
        for(ui i = 0; i < vertices_count_; i++) h[i] = deg[i] / graph.weight_[i], id[i] = 0;
        heap.init(h);
        double low = graph.subgraph_density_lower_bound;
        for(ui i = 0; i < vertices_count_; i++){
            HeapElement tmp = heap.pop();
            h[tmp.id] = tmp.val;
            id[tmp.id] = 1;
            low = std::max(low, tmp.val / 2);
            for(auto v: graph.getNeighbors(tmp.id)){
                if(id[v] != 0) continue;
                heap.dec(v, 1.0 / graph.weight_[v], tmp.val);
            }
        }
        
        ui new_vertices_count_ = 0;
        for(ui i = 0; i < vertices_count_; i++){
            if(h[i] >= low) id[i] = new_vertices_count_++;
            else id[i] = vertices_count_ + 1;
        }
        Graph copy(0, new_vertices_count_);
        copy.weight_.resize(new_vertices_count_);
        copy.map.reserve(new_vertices_count_);
        copy.subgraph_density_lower_bound = std::max(low, graph.subgraph_density_lower_bound);
        new_vertices_count_ = 0;
        for(ui i = 0; i < vertices_count_; i++){
            if(id[i] != vertices_count_ + 1){
                copy.map[new_vertices_count_] = graph.map[i];
                copy.weight_[new_vertices_count_++] = graph.weight_[i];
            }
        }
        for(ui i = 0; i < vertices_count_; i++){
            if(id[i] == vertices_count_ + 1) continue;
            for(auto v: graph.getNeighbors(i)){
                if(id[v] ==  vertices_count_ + 1|| id[i] < id[v]) continue;
                copy.addUndirectedEdge(id[i], id[v]);
            }
        }
        std::vector<Alpha> new_alpha;
        for(ui i = 0; i < graph.getEdgesCount(); i++){
            if(id[lp.alpha[i].id_first] != vertices_count_ + 1 && id[lp.alpha[i].id_second] != vertices_count_ + 1){
                lp.alpha[i].id_first = id[lp.alpha[i].id_first];
                lp.alpha[i].id_second= id[lp.alpha[i].id_second];
                new_alpha.push_back(lp.alpha[i]);
            }
        }
        lp.alpha = new_alpha;
        if(lp.type_){
            std::vector<Alpha> new_beta;
            for(ui i = 0; i < graph.getEdgesCount(); i++){
                if(id[lp.beta[i].id_first] != vertices_count_ + 1 && id[lp.beta[i].id_second] != vertices_count_ + 1){
                    lp.beta[i].id_first = id[lp.beta[i].id_first];
                    lp.beta[i].id_second= id[lp.beta[i].id_second];
                    new_beta.push_back(lp.beta[i]);
                }
            }
            lp.beta = new_beta;
        }
        lp.nodes_count_ = copy.getVerticesCount();
        lp.edges_count_ = copy.getEdgesCount();
        for(ui i = 0; i < lp.nodes_count_; i++) lp.r[0][i] = 0;
        for(ui i = 0; i < lp.edges_count_; i++){
            lp.r[0][lp.alpha[i].id_first] += lp.alpha[i].weight_first / copy.weight_[lp.alpha[i].id_first];
            lp.r[0][lp.alpha[i].id_second] += lp.alpha[i].weight_second / copy.weight_[lp.alpha[i].id_second];
        }
        lp.weight = copy.weight_;
        graph = copy;
        graph = copy;
    }
}
void Reduction::UndirectedStableReduction(Graph &graph, LinearProgramming &lp){
    std::vector<ui> y;
    std::vector<std::pair<double, VertexID>> tmp;
    ui n = graph.getVerticesCount();
    ui m = graph.getEdgesCount();
    y.resize(n);
    tmp.resize(n);
    for (ui u = 0; u < n; u++) {
        y[u] = 0;
        tmp[u] = std::make_pair(-lp.r[0][u], u);
    }
    for (ui i = 0; i < m; i++) {
        if (lp.r[0][lp.alpha[i].id_first] < lp.r[0][lp.alpha[i].id_second] ||
            (lp.r[0][lp.alpha[i].id_first] == lp.r[0][lp.alpha[i].id_second] &&
             lp.alpha[i].id_first > lp.alpha[i].id_second))
            y[lp.alpha[i].id_first] += 1.0 / lp.weight[lp.alpha[i].id_first];
        else
            y[lp.alpha[i].id_second] += 1.0 / lp.weight[lp.alpha[i].id_second];
    }
    sort(tmp.begin(), tmp.end());
    std::vector<double> val;
    std::vector<ui> sz;
    std::vector<ui> belong;
    belong.resize(n);
    sz.resize(n);
    val.resize(n);
    ui top = -1;
    double eps = 1e-6;
    for (ui u = 0; u < n; u++) {
        top++;
        val[top] = y[tmp[u].second];
        sz[top] = 1;
        while (top >= 1 && (val[top] + val[top - 1]) / (sz[top] + sz[top - 1]) > val[top - 1] / sz[top - 1]) {
            val[top - 1] += val[top];
            sz[top - 1] += sz[top];
            top--;
        }
    }
    ui pos = 0;
    for (ui u = 0; u < n; u++) {
        if (!sz[pos]) pos++;
        sz[pos]--;
        belong[tmp[u].second] = pos;
    }
    for (ui i = 0; i < m; i++) {
        if (belong[lp.alpha[i].id_first] == belong[lp.alpha[i].id_second]) continue;
        if (belong[lp.alpha[i].id_first] < belong[lp.alpha[i].id_second]) {
            lp.r[0][lp.alpha[i].id_second] += lp.alpha[i].weight_first / lp.weight[lp.alpha[i].id_first];
            lp.r[0][lp.alpha[i].id_first] -= lp.alpha[i].weight_first / lp.weight[lp.alpha[i].id_second];
            lp.alpha[i].weight_first = 0;
            lp.alpha[i].weight_second = 1;
        }
        if (belong[lp.alpha[i].id_first] > belong[lp.alpha[i].id_second]) {
            lp.r[0][lp.alpha[i].id_first] += lp.alpha[i].weight_second / lp.weight[lp.alpha[i].id_first];
            lp.r[0][lp.alpha[i].id_second] -= lp.alpha[i].weight_second / lp.weight[lp.alpha[i].id_second];
            lp.alpha[i].weight_second = 0;
            lp.alpha[i].weight_first = 1;
        }
    }
    double minn = 1e20;
    std::vector<ui> rev;
    rev.resize(n);
    pos = 0;
    for (ui u = 0; u < n; u++) {
        rev[tmp[u].second] = u + 1;
        if (belong[tmp[u].second] > pos) {
            pos++;
        }
        if (pos != 0 && lp.r[0][tmp[u].second] >= minn) {
            std::cout<<minn<<" haha "<<lp.r[0][tmp[u].second]<<std::endl;
            return;
        } else if (pos == 0)
            minn = std::min(minn, lp.r[0][tmp[u].second]);
    }
    std::cout<<"haha"<<std::endl;
    std::vector<ui> id;
    id.resize(n);
    for(ui u = 0; u < n; u++) id[u] = 0;
    ui new_vertices_count_ = 0;
    for(ui i = 0; i < n; i++){
        if(belong[tmp[i].second] == 0){
            id[tmp[i].second] = 1 + new_vertices_count_;
            new_vertices_count_++;
        }
    }
    Graph copy(0, new_vertices_count_);
    copy.weight_.resize(new_vertices_count_);
    for(ui i = 0; i < new_vertices_count_; i++) copy.weight_[i] = 1;
    for(ui i = 0; i < n; i++){
        if(id[i] < 1) continue;
        for(auto v: graph.getNeighbors(i)){
            if(id[v] < 1 || id[i] < id[v]) continue;
            copy.addUndirectedEdge(id[i] - 1, id[v] - 1);
        }
    }
    std::vector<Alpha> new_alpha;
    for(ui i = 0; i < graph.getEdgesCount(); i++){
        if(id[lp.alpha[i].id_first] >= 1 && id[lp.alpha[i].id_second] >= 1){
            lp.alpha[i].id_first = id[lp.alpha[i].id_first] - 1;
            lp.alpha[i].id_second= id[lp.alpha[i].id_second] - 1;
            new_alpha.push_back(lp.alpha[i]);
        }
    }
    lp.alpha = new_alpha;
    if(lp.type_){
        std::vector<Alpha> new_beta;
        for(ui i = 0; i < graph.getEdgesCount(); i++){
            if(id[lp.beta[i].id_first] >= 1 && id[lp.beta[i].id_second] >= 1){
                lp.beta[i].id_first = id[lp.beta[i].id_first] - 1;
                lp.beta[i].id_second= id[lp.beta[i].id_second] - 1;
                new_beta.push_back(lp.beta[i]);
            }
        }
        lp.beta = new_beta;
    }
    lp.nodes_count_ = copy.getVerticesCount();
    lp.edges_count_ = copy.getEdgesCount();
    for(ui i = 0; i < lp.nodes_count_; i++) lp.r[0][i] = 0;
    for(ui i = 0; i < lp.edges_count_; i++){
        lp.r[0][lp.alpha[i].id_first] += lp.alpha[i].weight_first;
        lp.r[0][lp.alpha[i].id_second] += lp.alpha[i].weight_second;
    }
    graph = copy;
}
