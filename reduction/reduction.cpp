//
// Created by yy on 11/23/23.
//
#include "reduction.h"

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
                                bool &is_init, bool is_dc, bool is_divide_by_number, bool is_exact) {
    if (!is_init) {
        is_init = true;
        xycore.xyCoreInitialization(graph);
    }
    ui x, y;
    if (!is_dc) {
        double ratio_sqrt = sqrt(ratios.first / ratios.second);
        x = std::max(static_cast<int>(ceil(l / 2 / ratio_sqrt)), 1);
        y = std::max(static_cast<int>(ceil(ratio_sqrt * l / 2)), 1);
    } else {
        double ratio;
        if (is_divide_by_number){
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
        double ratio_right_sqrt = sqrt(std::min(ratios.second, ratio * 2));
        double ratio_left_sqrt = sqrt(std::max(ratios.first, ratio / 2));
//        double ratio_right_sqrt = sqrt(ratios.second);
//        double ratio_left_sqrt = sqrt(ratios.first);
        x = std::max(static_cast<int>(ceil(l / 2 / ratio_right_sqrt)), 1);
        y = std::max(static_cast<int>(ceil(ratio_left_sqrt * l / 2)), 1);
    }

    xycore.generateXYCore(graph, x_y_core, x, y, is_exact);
}

void Reduction::stableSetReduction(Graph &graph, LinearProgramming &lp,
                                   std::vector<std::pair<VertexID, VertexID>> edges, bool stable_set_reduction) {
    if (!stable_set_reduction)
        return;
    Graph stable_subgraph = Graph(true, graph.getVerticesCount());
    std::vector<bool> selected[2];
    for (ui i = 0; i < 2; i++)
        selected[i].resize(graph.getVerticesCount(), false);
    for (auto &edge: edges) {
        stable_subgraph.addDirectedEdge(edge.first, edge.second);
        selected[0][edge.first] = true;
        selected[1][edge.second] = true;
    }
    ui i = 0, j = lp.edges_count_ - 1;
    while (i != j) {
        while (i < j && selected[0][lp.alpha[i].id_first] && selected[1][lp.alpha[i].id_second])
            i++;
        while (j > i && (!selected[0][lp.alpha[j].id_first] || !selected[1][lp.alpha[j].id_second]))
            j--;
//        printf("%d, %d\n", i, j);
        std::swap(lp.alpha[i], lp.alpha[j]);
    }
    lp.edges_count_ = edges.size();
    lp.alpha.resize(lp.edges_count_);
    lp.sort(stable_subgraph);
//    for (ui i = 0; i < lp.edges_count_; i++) {
//        if (!lp.alpha[i].is_selected)
//            continue;
//        if (!selected[0][lp.alpha[i].id_first] || !selected[1][lp.alpha[i].id_second]) {
//            lp.alpha[i].is_selected = false;
//        }
//    }
    graph = stable_subgraph;
}

void Reduction::wCoreReduction(Graph &graph, WCore &w_core) {
    w_core.generateMaxWCore(graph);
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
//VertexID i = 0, j = 0;
//while(deg[0][vert_[0][i]] < x || deg[1][vert_[1][j]] < y) {
//for (; deg[0][vert_[0][i]] < x; i++) {
//auto u = vert_[0][i];
//for (auto v : graph.getOutNeighbors(u)) {
//if (deg[1][v] >= y) {
////                    dec_deg(1, v);
//auto dt = deg[1][v], pt = pos_[1][v];
//auto pw = bin_[1][dt], w = vert_[1][pw];
//if (v != w) {
//pos_[1][v] = pw; vert_[1][pt] = w;
//pos_[1][w] = pt; vert_[1][pw] = v;
//}
//++bin_[1][dt];
//--deg[1][v];
//}
//}
//}
//
//for (; deg[1][vert_[1][j]] < y; j++) {
//auto u = vert_[1][j];
//for (auto v : graph.getInNeighbors(u)) {
//if (deg[0][v] >= x) {
////                    dec_deg(0, v);
//auto dt = deg[0][v], pt = pos_[0][v];
//auto pw = bin_[0][dt], w = vert_[0][pw];
//if (v != w) {
//pos_[0][v] = pw; vert_[0][pt] = w;
//pos_[0][w] = pt; vert_[0][pw] = v;
//}
//++bin_[0][dt];
//--deg[0][v];
//}
//}
//}
//}
//auto *vertices = new std::vector<VertexID>[2];
//std::vector<std::pair<VertexID, VertexID>> edges;
//for (; i < vertices_count; i++) {
//auto u = vert_[0][i];
//vertices[0].push_back(u);
//for (auto v : graph.getOutNeighbors(u)) {
//if (deg[1][v] >= y) {
//vertices[1].push_back(v);
//}
//}
//}
//sort( vertices[1].begin(), vertices[1].end() );
//vertices[1].erase( unique( vertices[1].begin(), vertices[1].end() ), vertices[1].end() );
//delete[] deg;
//return vertices;