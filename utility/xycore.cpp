//
// Created by yy on 12/16/23.
//

#include "xycore.h"

XYCore::XYCore() {}

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

void XYCore::generateXYCore(const Graph &graph, Graph &x_y_core, ui x, ui y) {
    auto n = graph.getVerticesCount();

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

inline void XYCore::decDeg(ui cur, VertexID t) {
    ui dt = degrees[cur][t], pt = pos[cur][t];
    ui pw = bin[cur][dt];
    VertexID w = vert[cur][pw];
    if (t != w) {
        pos[cur][t] = pw; vert[cur][pt] = w;
        pos[cur][w] = pt; vert[cur][pw] = t;
    }
    ++bin[cur][dt];
    --degrees[cur][t];
}

ui XYCore::getDelta(const Graph &graph) {
    ui vertices_count = graph.getVerticesCount();
    ui i = bin[0][1], j = bin[1][1];
    ui delta = 0;

    while(i < vertices_count && j < vertices_count){
        ui cur;
        VertexID u;

        if(degrees[0][vert[0][i]] <= degrees[1][vert[1][j]]) {
            cur = 0;
            u = vert[0][i];
            i++;
        } else {
            cur = 1;
            u = vert[1][j];
            j++;
        }
        delta = std::max(delta, degrees[cur][u]);
        for (auto v : (cur) ? graph.getInNeighbors(u) : graph.getOutNeighbors(u)) {
            if (degrees[1 - cur][v] > degrees[cur][u]) {
                decDeg(1 - cur, v);
            }
        }
    }

    return delta;
}

ui XYCore::skyline_core_num(Graph &graph, ui cur, ui x, ui y, bool reduced) {
    ui vertices_count = graph.getVerticesCount();

    if (y > std::upper_bound(bin[1 - cur].begin(), bin[1 - cur].end(), vertices_count - x) - bin[1 - cur].begin() - 1)
        return 0;

    if(reduced) {
        std::vector<bool> flag(static_cast<unsigned long>(vertices_count), false);
        ui ret = 0;
        for (int i = 0; i < vertices_count; i++) {
            VertexID u = vert[1 - cur][i];
            ret = std::max(ret, degrees[1 - cur][u]);
            for (auto v : (1 - cur) ? graph.getInNeighbors(u) : graph.getOutNeighbors(u)) {
                if (flag[v]) continue;
                if (--degrees[cur][v] < x) {
                    flag[v] = true;
                    for (auto t : (cur) ? graph.getInNeighbors(v) : graph.getOutNeighbors(v)) {
                        if (t != u && degrees[1 - cur][t] > degrees[1 - cur][u]) {
                            decDeg(1 - cur, t);
                        }
                    }
                }
            }
        }
        return ret;
    } else{
        Graph reduced_graph(vertices_count);
        generateXYCore(graph, reduced_graph, x, y);
        auto ret = skyline_core_num(reduced_graph, 0, x, y, true);
        return ret;
    }
}