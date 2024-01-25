//
// Created by yy on 12/16/23.
//

#include "xycore.h"

XYCore::XYCore() {
    vert.resize(2);
    bin.resize(2);
    pos.resize(2);
    degrees.resize(2);
}

void XYCore::xyCoreInitialization(const Graph &graph) {


    auto vertices_count = graph.getVerticesCount();
    degrees[0] = graph.getOutDegrees();
    degrees[1] = graph.getInDegrees();
    for (int i = 0; i < 2; i++) {
        vert[i].clear();
        bin[i].clear();
        pos[i].clear();

        vert[i].resize(vertices_count, 0);
        bin[i].resize(vertices_count + 1, 0);
        pos[i].resize(vertices_count, 0);
        for (VertexID v = 0; v < vertices_count; v++) {
            ++bin[i][degrees[i][v]];
        }
        ui start = 0;
        for (int d = 0; d <= vertices_count; d++) {
            ui num = bin[i][d];
            bin[i][d] = start;
            start += num;
        }
        for (VertexID v = 0; v < vertices_count; v++) {
            pos[i][v] = bin[i][degrees[i][v]]++;
            vert[i][pos[i][v]] = v;
        }
        for (int d = vertices_count; d > 0; d--) {
            bin[i][d] = bin[i][d - 1];
        }
        bin[i][0] = 0;
    }
}

void XYCore::generateXYCore(const Graph &graph, Graph &x_y_core, ui x, ui y, bool is_exact) {
    auto n = graph.getVerticesCount();
    if (is_exact) {
        std::vector<std::vector<VertexID>> vert_copy(2);
        std::vector<std::vector<ui>> bin_copy(2);
        std::vector<std::vector<ui>> pos_copy(2);
        std::vector<std::vector<ui>> degrees_copy(2);

        for (ui i = 0; i < 2; i++) {
            vert_copy[i] = vert[i];
            bin_copy[i] = bin[i];
            pos_copy[i] = pos[i];
            degrees_copy[i] = degrees[i];
        }
        VertexID i = 0, j = 0;
        while (degrees_copy[0][vert_copy[0][i]] < x || degrees_copy[1][vert_copy[1][j]] < y) {
            for (; i < n && degrees_copy[0][vert_copy[0][i]] < x; i++) {
                auto u = vert_copy[0][i];
                for (auto v: graph.getOutNeighbors(u)) {
                    if (degrees_copy[1][v] >= y) {
                        decDeg(1, v, vert_copy, bin_copy, pos_copy, degrees_copy);
                    }
                }
            }

            for (; j < n && degrees_copy[1][vert_copy[1][j]] < y; j++) {
                auto u = vert_copy[1][j];
                for (auto v: graph.getInNeighbors(u)) {
                    if (degrees_copy[0][v] >= x) {
                        decDeg(0, v, vert_copy, bin_copy, pos_copy, degrees_copy);
                    }
                }
            }

            if (i == n || j == n)
                break;
        }
//    std::vector<bool> is_selected[2];
//    for (ui i = 0; i < 2; i++)
//        is_selected[i].resize(n, false);
//    std::vector<std::pair<VertexID, VertexID>> edges;
        for (; i < n; i++) {
            auto u = vert_copy[0][i];
//        vertices[0].push_back(u);
            for (auto v: graph.getOutNeighbors(u)) {
                if (degrees_copy[1][v] >= y) {
//                vertices[1].push_back(v);
                    x_y_core.addDirectedEdge(u, v);
                }
            }
        }
    } else {
        int cur = 0;
        if (x > *std::max_element(degrees[cur].begin(), degrees[cur].end())) {
            return;
        }
        for (ui i = bin[cur][x]; i < n; i++) {
            VertexID u = vert[cur][i];
            auto adj = graph.getOutNeighbors(u);
            std::sort(adj.begin(), adj.end(),
                      [&](const int &a, const int &b) -> bool {
                          return degrees[1 - cur][a] > degrees[1 - cur][b];
                      });
            VertexID v = adj[std::max(int(x) - 1, 0)];
            if (degrees[1 - cur][v] < y) continue;
            for (auto t: adj) {
                if (degrees[1 - cur][t] < y)
                    break;
                x_y_core.addDirectedEdge(u, t);
            }
        }
    }
//    sort( vertices[1].begin(), vertices[1].end() );
//    vertices[1].erase( unique( vertices[1].begin(), vertices[1].end() ), vertices[1].end() );
//    delete[] degrees_;
//    return vertices;
}

inline void XYCore::decDeg(ui cur, VertexID t, std::vector<std::vector<VertexID>> &vert_copy,
                           std::vector<std::vector<ui>> &bin_copy, std::vector<std::vector<ui>> &pos_copy,
                           std::vector<std::vector<ui>> &degrees_copy) {
    ui dt = degrees_copy[cur][t], pt = pos_copy[cur][t];
    ui pw = bin_copy[cur][dt];
    VertexID w = vert_copy[cur][pw];
    if (t != w) {
        pos_copy[cur][t] = pw;
        vert_copy[cur][pt] = w;
        pos_copy[cur][w] = pt;
        vert_copy[cur][pw] = t;
    }
    ++bin_copy[cur][dt];
    --degrees_copy[cur][t];
//    int dt = deg[cur][t], pt = pos[cur][t];
//    int pw = bin[cur][dt], w = vert[cur][pw];
//    if (t != w) {
//        pos[cur][t] = pw; vert[cur][pt] = w;
//        pos[cur][w] = pt; vert[cur][pw] = t;
//    }
//    ++bin[cur][dt];
//    --deg[cur][t];
}

ui XYCore::getDelta(const Graph &graph) {
    ui vertices_count = graph.getVerticesCount();
    std::vector<std::vector<VertexID>> vert_copy(2);
    std::vector<std::vector<ui>> bin_copy(2);
    std::vector<std::vector<ui>> pos_copy(2);
    std::vector<std::vector<ui>> degrees_copy(2);

    for (ui i = 0; i < 2; i++) {
        vert_copy[i] = vert[i];
        bin_copy[i] = bin[i];
        pos_copy[i] = pos[i];
        degrees_copy[i] = degrees[i];
    }
    ui i = bin_copy[0][1], j = bin_copy[1][1];
    ui delta = 0;

    while (i < vertices_count && j < vertices_count) {
        ui cur;
        VertexID u;

        if (degrees_copy[0][vert_copy[0][i]] <= degrees_copy[1][vert_copy[1][j]]) {
            cur = 0;
            u = vert_copy[0][i];
            i++;
        } else {
            cur = 1;
            u = vert_copy[1][j];
            j++;
        }
        delta = std::max(delta, degrees_copy[cur][u]);
        for (auto v: (cur) ? graph.getInNeighbors(u) : graph.getOutNeighbors(u)) {
            if (degrees_copy[1 - cur][v] > degrees_copy[cur][u]) {
                decDeg(1 - cur, v, vert_copy, bin_copy, pos_copy, degrees_copy);
            }
        }
    }

    return delta;
}

ui XYCore::skyline_core_num(Graph &graph, ui cur, ui x, ui y, bool reduced) {
    ui vertices_count = graph.getVerticesCount();
    std::vector<std::vector<VertexID>> vert_copy(2);
    std::vector<std::vector<ui>> bin_copy(2);
    std::vector<std::vector<ui>> pos_copy(2);
    std::vector<std::vector<ui>> degrees_copy(2);

    for (ui i = 0; i < 2; i++) {
        vert_copy[i] = vert[i];
        bin_copy[i] = bin[i];
        pos_copy[i] = pos[i];
        degrees_copy[i] = degrees[i];
    }
    if (y > std::upper_bound(bin_copy[1 - cur].begin(), bin_copy[1 - cur].end(), vertices_count - x) -
            bin_copy[1 - cur].begin() - 1)
        return 0;

    if (reduced) {
        std::vector<bool> flag(static_cast<unsigned long>(vertices_count), false);
        ui ret = 0;
        for (int i = 0; i < vertices_count; i++) {
            VertexID u = vert_copy[1 - cur][i];
            ret = std::max(ret, degrees_copy[1 - cur][u]);
            for (auto v: (1 - cur) ? graph.getInNeighbors(u) : graph.getOutNeighbors(u)) {
                if (flag[v]) continue;
                if (--degrees_copy[cur][v] < x) {
                    flag[v] = true;
                    for (auto t: (cur) ? graph.getInNeighbors(v) : graph.getOutNeighbors(v)) {
                        if (t != u && degrees_copy[1 - cur][t] > degrees_copy[1 - cur][u]) {
                            decDeg(1 - cur, t, vert_copy, bin_copy, pos_copy, degrees_copy);
                        }
                    }
                }
            }
        }
        return ret;
    } else {
        Graph reduced_graph(true, vertices_count);
        generateXYCore(graph, reduced_graph, x, y);
        XYCore core;
        core.xyCoreInitialization(reduced_graph);
        auto ret = core.skyline_core_num(reduced_graph, cur, x, y, true);
        return ret;
    }
}