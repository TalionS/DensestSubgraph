//
// Created by yy on 12/16/23.
//

#include "xycore.h"

XYCore::XYCore() {
    vert.resize(2);
    bin.resize(2);
    pos.resize(2);
    degrees.resize(2);
    adj.resize(2);
    max_degrees.resize(2);
};

XYCore &XYCore::operator=(const XYCore &other) {
    if (this != &other) {
        vert = other.vert;
        bin = other.bin;
        pos = other.pos;
        degrees = other.degrees;
        adj = other.adj;
    }
    return *this;
};

void XYCore::xyCoreInitialization(Graph &graph, bool sort) {


    auto vertices_count = graph.getVerticesCount();
    degrees[0] = graph.getOutDegrees();
    degrees[1] = graph.getInDegrees();
    adj = graph.getAdjList();
    if (sort)
        for (ui u = 0; u < vertices_count; u++) {
            std::sort(adj[0][u].begin(), adj[0][u].end(),
                      [&](const int &a, const int &b) -> bool {
                          return degrees[1][a] > degrees[1][b];
                      });
        }
    for (int i = 0; i < 2; i++) {
        max_degrees[i] = *std::max_element(degrees[i].begin(), degrees[i].end());
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


void XYCore::generateXYCore(const Graph &graph, Graph &subgraph, ui x, ui y, bool is_exact, bool is_map, bool is_copy) {
    auto n = graph.getVerticesCount();
    if (is_exact) {
        VertexID i = 0, j = 0;
        while (degrees[0][vert[0][i]] < x || degrees[1][vert[1][j]] < y) {
            for (; i < n && degrees[0][vert[0][i]] < x; i++) {
                auto u = vert[0][i];
                for (auto v: adj[0][u]) {
                    if (degrees[1][v] >= y) {
                        decDeg(1, v, vert, bin, pos, degrees);
                    }
                }
            }

            for (; j < n && degrees[1][vert[1][j]] < y; j++) {
                auto u = vert[1][j];
                for (auto &v: adj[1][u]) {
                    if (degrees[0][v] >= x) {
                        decDeg(0, v, vert, bin, pos, degrees);
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
        if (!is_map) {
            subgraph.setVerticesCount(n);
            for (; i < n; i++) {
                auto u = vert[0][i];
//        vertices[0].push_back(u);
                for (auto &v: adj[0][u]) {
                    if (degrees[1][v] >= y) {
//                vertices[1].push_back(v);
                        subgraph.addDirectedEdge(u, v);
                    }
                }
            }
        } else {
            std::vector<bool> is_selected(n, false);
            std::vector<std::pair<VertexID, VertexID>> edges;
            ui vertex_id = 0;
            for (; i < n; i++) {
                auto u = vert[0][i];
                is_selected[u] = true;
                for (auto &v: adj[0][u]) {
                    edges.emplace_back(u, v);
                    is_selected[v] = true;
                }
            }
            std::vector<VertexID> map(n, 0);
            std::vector<VertexID> reverse_map;
            if (graph.is_mapped) {
                for (ui u = 0; u < n; u++)
                    if (is_selected[u]) {
                        reverse_map.emplace_back(graph.map[u]);
                        map[u] = vertex_id++;
                    }
            } else {
                for (ui u = 0; u < n; u++)
                    if (is_selected[u]) {
                        reverse_map.emplace_back(u);
                        map[u] = vertex_id++;
                    }
            }
            subgraph.setVerticesCount(vertex_id);
            for (auto &edge: edges)
                subgraph.addDirectedEdge(map[edge.first], map[edge.second]);
            subgraph.is_mapped = true;
            subgraph.map = reverse_map;
        }
        return;
    } else if (!is_map) {
        if (is_copy) {
            Graph copy(true, n);
            int cur = 0;
            if (x > max_degrees[cur]) {
                return;
            }
            for (ui i = bin[cur][x]; i < n; i++) {
                VertexID u = vert[cur][i];
                VertexID v = adj[0][u][std::max(int(x) - 1, 0)];
                if (degrees[1 - cur][v] < y) continue;
                for (auto &t: adj[0][u]) {
                    if (degrees[1 - cur][t] < y)
                        break;
                    copy.addDirectedEdge(u, t);
                }
            }
            if (graph.is_mapped) {
                copy.is_mapped = true;
                copy.map = graph.map;
            }
            subgraph = copy;
//        std::vector<std::vector<ui>> deg(2);
//        deg[0] = subgraph.getOutDegrees();
//        deg[1] = subgraph.getInDegrees();
//        for (VertexID u = 0; u < n; u++) {
//            if (deg[0][u])
//                subgraph.vertex_ids[0].emplace_back(u);
//            if (deg[1][u])
//                subgraph.vertex_ids[1].emplace_back(u);
//        }
            return;
        } else {
            subgraph.setVerticesCount(n);
            int cur = 0;
            if (x > max_degrees[cur]) {
                return;
            }
            for (ui i = bin[cur][x]; i < n; i++) {
                VertexID u = vert[cur][i];
                VertexID v = adj[0][u][std::max(int(x) - 1, 0)];
                if (degrees[1 - cur][v] < y) continue;
                for (auto &t: adj[0][u]) {
                    if (degrees[1 - cur][t] < y)
                        break;
                    subgraph.addDirectedEdge(u, t);
                }
            }
//        std::vector<std::vector<ui>> deg(2);
//        deg[0] = subgraph.getOutDegrees();
//        deg[1] = subgraph.getInDegrees();
//        for (VertexID u = 0; u < n; u++) {
//            if (deg[0][u])
//                subgraph.vertex_ids[0].emplace_back(u);
//            if (deg[1][u])
//                subgraph.vertex_ids[1].emplace_back(u);
//        }
            return;
        }
    } else {
        if (is_copy) {
            Graph copy(true, 0);
            int cur = 0;
            if (x > max_degrees[cur]) {
//        if (x  > *std::max_element(degrees[cur].begin(), degrees[cur].end())) {
                return;
            }
            std::vector<bool> is_selected(n, false);
            std::vector<std::pair<VertexID, VertexID>> edges;
            for (ui i = bin[cur][x]; i < n; i++) {
                VertexID u = vert[cur][i];
//            auto adj = graph.getOutNeighbors(u);
//            std::sort(adj.begin(), adj.end(),
//                      [&](const int &a, const int &b) -> bool {
//                          return degrees[1 - cur][a] > degrees[1 - cur][b];
//                      });
                VertexID v = adj[0][u][std::max(int(x) - 1, 0)];
                if (degrees[1 - cur][v] < y) continue;
                is_selected[u] = true;
                for (auto &t: adj[0][u]) {
                    if (degrees[1 - cur][t] < y)
                        break;
                    is_selected[t] = true;
                    edges.emplace_back(u, t);
//                x_y_core.addDirectedEdge(u, t);
                }
            }
            ui vertex_id = 0;
            std::vector<VertexID> map(n);
            std::vector<VertexID> reverse_map;
            if (graph.is_mapped) {
                for (ui u = 0; u < n; u++)
                    if (is_selected[u]) {
                        reverse_map.emplace_back(graph.map[u]);
                        map[u] = vertex_id++;
                    }
            } else {
                for (ui u = 0; u < n; u++)
                    if (is_selected[u]) {
                        reverse_map.emplace_back(u);
                        map[u] = vertex_id++;
                    }
            }
            copy.setVerticesCount(vertex_id);
            copy.is_mapped = true;
            copy.map = reverse_map;
            for (auto &edge: edges)
                copy.addDirectedEdge(map[edge.first], map[edge.second]);
            subgraph = copy;
            return;
        } else {
            int cur = 0;
            if (x > max_degrees[cur]) {
//        if (x  > *std::max_element(degrees[cur].begin(), degrees[cur].end())) {
                return;
            }
            std::vector<bool> is_selected(n, false);
            std::vector<std::pair<VertexID, VertexID>> edges;
            for (ui i = bin[cur][x]; i < n; i++) {
                VertexID u = vert[cur][i];
//            auto adj = graph.getOutNeighbors(u);
//            std::sort(adj.begin(), adj.end(),
//                      [&](const int &a, const int &b) -> bool {
//                          return degrees[1 - cur][a] > degrees[1 - cur][b];
//                      });
                VertexID v = adj[0][u][std::max(int(x) - 1, 0)];
                if (degrees[1 - cur][v] < y) continue;
                is_selected[u] = true;
                for (auto &t: adj[0][u]) {
                    if (degrees[1 - cur][t] < y)
                        break;
                    is_selected[t] = true;
                    edges.emplace_back(u, t);
//                x_y_core.addDirectedEdge(u, t);
                }
            }
            ui vertex_id = 0;
            std::vector<VertexID> map(n);
            std::vector<VertexID> reverse_map;
            if (graph.is_mapped) {
                for (ui u = 0; u < n; u++)
                    if (is_selected[u]) {
                        reverse_map.emplace_back(graph.map[u]);
                        map[u] = vertex_id++;
                    }
            } else {
                for (ui u = 0; u < n; u++)
                    if (is_selected[u]) {
                        reverse_map.emplace_back(u);
                        map[u] = vertex_id++;
                    }
            }
            subgraph.setVerticesCount(vertex_id);
            subgraph.is_mapped = true;
            subgraph.map = reverse_map;
            for (auto &edge: edges)
                subgraph.addDirectedEdge(map[edge.first], map[edge.second]);
//        subgraph = x_y_core;
            return;
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
    std::vector<std::vector<VertexID>> vert_copy = vert;
    std::vector<std::vector<ui>> bin_copy = bin;
    std::vector<std::vector<ui>> pos_copy = pos;
    std::vector<std::vector<ui>> degrees_copy = degrees;

//    for (ui i = 0; i < 2; i++) {
//        vert_copy[i] = vert[i];
//        bin_copy[i] = bin[i];
//        pos_copy[i] = pos[i];
//        degrees_copy[i] = degrees[i];
//    }
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
        for (auto &v: adj[cur][u]) {
            if (degrees_copy[1 - cur][v] > degrees_copy[cur][u]) {
                decDeg(1 - cur, v, vert_copy, bin_copy, pos_copy, degrees_copy);
            }
        }
    }

    return delta;
}

ui XYCore::skyline_core_num(Graph &graph, ui cur, ui x, ui y, bool reduced) {
    ui vertices_count = graph.getVerticesCount();
//    std::vector<std::vector<VertexID>> vert(2);
//    std::vector<std::vector<ui>> bin(2);
//    std::vector<std::vector<ui>> pos(2);
//    std::vector<std::vector<ui>> degrees(2);
//
//    for (ui i = 0; i < 2; i++) {
//        vert[i] = vert[i];
//        bin[i] = bin[i];
//        pos[i] = pos[i];
//        degrees[i] = degrees[i];
//    }
    if (y > std::upper_bound(bin[1 - cur].begin(), bin[1 - cur].end(), vertices_count - x) -
            bin[1 - cur].begin() - 1)
        return 0;

    if (reduced) {
        std::vector<bool> flag(static_cast<unsigned long>(vertices_count), false);
        ui ret = 0;
//        auto adj = graph.getAdjList();
        for (ui i = bin[1 - cur][1]; i < vertices_count; i++) {
            VertexID u = vert[1 - cur][i];
            ret = std::max(ret, degrees[1 - cur][u]);
            for (auto &v: adj[1 - cur][u]) {
                if (flag[v]) continue;
                if (--degrees[cur][v] < x) {
                    flag[v] = true;
                    for (auto &t: adj[cur][v]) {
                        if (t != u && degrees[1 - cur][t] > degrees[1 - cur][u]) {
                            decDeg(1 - cur, t, vert, bin, pos, degrees);
                        }
                    }
                }
            }
        }
        return ret;
    } else {
//        auto begin = clock();
        Graph reduced_graph(true, 0);
//        auto end = clock();
//        printf("construction time: %f\n", (double) (end - begin) / CLOCKS_PER_SEC);
        generateXYCore(graph, reduced_graph, x, y, false, true, false);
//        auto reduced_graph = generateXYCore(graph, x, y, false, true);
//        printf("%d\n", reduced_graph.getVerticesCount());
        if (!reduced_graph.getVerticesCount())
            return 0;
        XYCore core;
        core.xyCoreInitialization(reduced_graph);
//        auto end = clock();
//        printf("construction time: %f\n", (double) (end - begin) / CLOCKS_PER_SEC);
        auto ret = core.skyline_core_num(reduced_graph, cur, x, y, true);
//        auto end2 = clock();
//        printf("calculation time: %f\n", (double) (end2 - end) / CLOCKS_PER_SEC);
        return ret;
    }
}