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
//        assert(degree_start[degrees[u]] == i);
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

void Reduction::xyCoreReduction(Graph &graph, Graph &x_y_core, std::pair<double, double> ratio, double &l, double &r, bool &is_init, bool is_dc) {
    if(!is_init) {
        is_init = true;
        xycore.xyCoreInitialization(graph);
    }
    ui x, y;
    if(!is_dc) {
        double ratio_sqrt = sqrt(ratio.first / ratio.second);
        x = std::max(static_cast<int>(ceil(l / 2 / ratio_sqrt)), 1);
        y = std::max(static_cast<int>(ceil(ratio_sqrt * l / 2)), 1);
    }
    else {
        double ratio_right_sqrt = sqrt(ratio.second);
        double ratio_left_sqrt = sqrt(ratio.first);
        x = std::max(static_cast<int>(ceil(l / 2 / ratio_right_sqrt)), 1);
        y = std::max(static_cast<int>(ceil(ratio_left_sqrt * l / 2)), 1);
    }

    xycore.generateXYCore(graph, x_y_core, x, y);
}

void Reduction::kCoreReduction(Graph &graph, double &l, double &r){
    ui vertices_count = graph.getVerticesCount();
    std::vector<ui> degrees = graph.getDegrees();
    std::vector<ui> newid(vertices_count, 0);
    std::queue<VertexID> q;
    for(int i = 0; i < vertices_count; i++){
        if(degrees[i] < l) q.push(i), newid[i] = vertices_count;
    }
    while(q.size()){
        int u = q.front();
        q.pop();
        for (ui v: graph.getNeighbors(u)) {
            if(newid[v] == vertices_count) continue;
            if ((--degrees[v]) < l) {
                q.push(v);
                newid[v] = vertices_count;
            }
        }
    }
    ui new_vertices_count = 0;
    for(int i = 0; i < vertices_count; i++){
        if(newid[i] != vertices_count) newid[i] = new_vertices_count++;
    }
    Graph kcore = Graph(false, new_vertices_count);
    for(int i = 0; i < vertices_count; i++){
        if(newid[i] == vertices_count) continue;
        for (ui j: graph.getNeighbors(i)) {
            if(newid[j] == vertices_count) continue;
            if(newid[i] < newid[j]) kcore.addUndirectedEdge(newid[i], newid[j]);
        }
    }
    kcore.subgraph_density_lower_bound = l;
    kcore.subgraph_density_upper_bound = r;
    graph = kcore;
    auto dgrees = graph.getDegrees();
}
//VertexID i = 0, j = 0;
//while(degrees[0][vert[0][i]] < x || degrees[1][vert[1][j]] < y) {
//for (; degrees[0][vert[0][i]] < x; i++) {
//auto u = vert[0][i];
//for (auto v : graph.getOutNeighbors(u)) {
//if (degrees[1][v] >= y) {
////                    dec_deg(1, v);
//auto dt = degrees[1][v], pt = pos[1][v];
//auto pw = bin[1][dt], w = vert[1][pw];
//if (v != w) {
//pos[1][v] = pw; vert[1][pt] = w;
//pos[1][w] = pt; vert[1][pw] = v;
//}
//++bin[1][dt];
//--degrees[1][v];
//}
//}
//}
//
//for (; degrees[1][vert[1][j]] < y; j++) {
//auto u = vert[1][j];
//for (auto v : graph.getInNeighbors(u)) {
//if (degrees[0][v] >= x) {
////                    dec_deg(0, v);
//auto dt = degrees[0][v], pt = pos[0][v];
//auto pw = bin[0][dt], w = vert[0][pw];
//if (v != w) {
//pos[0][v] = pw; vert[0][pt] = w;
//pos[0][w] = pt; vert[0][pw] = v;
//}
//++bin[0][dt];
//--degrees[0][v];
//}
//}
//}
//}
//auto *vertices = new std::vector<VertexID>[2];
//std::vector<std::pair<VertexID, VertexID>> edges;
//for (; i < vertices_count; i++) {
//auto u = vert[0][i];
//vertices[0].push_back(u);
//for (auto v : graph.getOutNeighbors(u)) {
//if (degrees[1][v] >= y) {
//vertices[1].push_back(v);
//}
//}
//}
//sort( vertices[1].begin(), vertices[1].end() );
//vertices[1].erase( unique( vertices[1].begin(), vertices[1].end() ), vertices[1].end() );
//delete[] degrees;
//return vertices;