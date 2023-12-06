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
};

void Reduction::generateXYCore(const Graph &graph, Graph &x_y_core, ui x, ui y) {
    auto n = graph.getVerticesCount();
    auto *degrees = new std::vector<ui>[2];
    degrees[0] = graph.getOutDegrees();
    degrees[1] = graph.getInDegrees();

    int cur = 0;
    if (x > *std::max_element(degrees[cur].begin(), degrees[cur].end())) {
        return;
    }
    for (ui i = bin[cur][x]; i < n; i++) {
        VertexID u = vert[cur][i];
        auto adj = graph.getOutNeighbors(u);
        std::sort(adj.begin(), adj.end(),
                  [&](const int& a, const int& b) -> bool
                  {
                      return degrees[1 - cur][a] > degrees[1 - cur][b];
                  });
        VertexID v = adj[std::max(int(x) - 1, 0)];
        if (degrees[1 - cur][v] < y) continue;
        for (auto t :adj) {
            if (degrees[1 - cur][t] < y)
                break;
            x_y_core.addDirectedEdge(u, t);
        }
    }

}

void Reduction::xyCoreInitialization(const Graph &graph) {
    auto vertices_count = graph.getVerticesCount();
    auto *degrees = new std::vector<ui>[2];
    degrees[0] = graph.getOutDegrees();
    degrees[1] = graph.getInDegrees();
    vert = new std::vector<VertexID>[2];
    bin = new std::vector<ui>[2];
    pos = new std::vector<ui>[2];
    for (int i = 0; i < 2; i++) {
        bin[i].resize(vertices_count + 1, 0);
        vert[i].resize(vertices_count, 0);
        pos[i].resize(vertices_count, 0);
        for(VertexID v = 0; v < vertices_count; v++){
            ++bin[i][degrees[i][v]];
        }
        ui start = 0;
        for(int d = 0; d <= vertices_count; d++){
            ui num = bin[i][d];
            bin[i][d] = start;
            start += num;
        }
        for(VertexID v = 0; v < vertices_count; v++){
            pos[i][v] = bin[i][degrees[i][v]]++;
            vert[i][pos[i][v]] = v;
        }
        for(int d = vertices_count; d > 0; d--){
            bin[i][d] = bin[i][d - 1];
        }
        bin[i][0] = 0;
    }
}

void Reduction::xyCoreReduction(Graph &graph, Graph &x_y_core, std::pair<double, double> ratio, double &l, double &r, bool &is_init, bool is_dc) {
    if(!is_init) {
        is_init = true;
        xyCoreInitialization(graph);
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

    generateXYCore(graph, x_y_core, x, y);
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