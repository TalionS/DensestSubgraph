#include "wcore.h"

void WCore::generateMaxWCore(Graph &graph, Graph &subgraph) {
    auto n = graph.getVerticesCount();
    auto adj = graph.getAdjList();
    std::vector<std::vector<ui>> deg(2);
//    std::vector<VertexID> vertices[2];
    deg[0] = graph.getOutDegrees();
    deg[1] = graph.getInDegrees();
    w = std::max(*std::max_element(deg[0].begin(), deg[0].end()),
                 *std::max_element(deg[1].begin(), deg[1].end())) - 1;
//    w = 27;

    for (ui i = 0; i < 2; i++)

        for (VertexID u = 0; u < n; u++)
            if (deg[i][u])
                vertices[i].push_back(u);
//    printf("S: %d, T: %d\n", vertices[0].size(), vertices[1].size());

    bool changed;
    while (!vertices[0].empty() && !vertices[1].empty()) {
        changed = true;
//        std::vector<std::pair<VertexID, VertexID>> peeled_edges;
        while (changed) {
            changed = false;
//            auto begin = std::chrono::steady_clock::now();
            for (auto u: vertices[0]) {
                ui d = deg[0][u];
                ui dout = 0;
//                auto out = graph.getOutNeighbors(u);
                std::vector<VertexID> new_neighbors;
                for (auto v: adj[0][u]) {
                    if (d * deg[1][v] > w) {
                        dout++;
                        new_neighbors.push_back(v);
                    } else {
//                        peeled_edges.emplace_back(u, v);
                    }
                }
                if (dout < deg[0][u]) {
                    deg[0][u] = dout;
                    adj[0][u] = new_neighbors;
                    changed = true;
                }
            }
//            auto end = std::chrono::steady_clock::now();
//            printf("%f\n", std::chrono::duration<double>(end - begin).count());

            for (auto v: vertices[1]) {
                ui d = deg[1][v];
                ui din = 0;
                std::vector<VertexID> new_neighbors;
                for (auto u: adj[1][v])
                    if (d * deg[0][u] > w) {
                        din++;
                        new_neighbors.push_back(u);
                    }
                if (din < deg[1][v]) {
                    deg[1][v] = din;
                    adj[1][v] = new_neighbors;
                    changed = true;
                }
            }
        }
//        for (auto edge: peeled_edges) {
//            edges.emplace_back(edge.first, edge.second);
//            induce_numbers.emplace_back(w);
//        }
        std::vector<ui> cnt(2, 0);
        for (ui i = 0; i < 2; i++)
            for (VertexID u = 0; u < n; u++)
                if (deg[i][u])
                    cnt[i]++;

        if (!cnt[0] || !cnt[1])
            break;

//        for (VertexID u = 0; u < n; u++)
//            if (deg[i][])

//        for (ui i = 0; i < 2; i++)
//            for(VertexID u = 0; u < n; u++)
//                if (!deg[i][u] && degrees[i][u])
//                    w_max.push_back(std::make_pair(w, ))

        degrees[0] = deg[0];
        degrees[1] = deg[1];
        vertices[0].clear();
        vertices[1].clear();

        for (ui i = 0; i < 2; i++)
            for (VertexID u = 0; u < n; u++)
                if (deg[i][u])
                    vertices[i].push_back(u);
//        printf("S: %d, T: %d\n", vertices[0].size(), vertices[1].size());

        ui min = UINT_MAX;
        for (auto u: vertices[0]) {
            ui d = deg[0][u];
            for (auto v: adj[0][u])
                if (min > d * deg[1][v])
                    min = d * deg[1][v];
        }

//        ui cur = vertices[0].size() > vertices[1].size() ? 1 : 0;
//        ui max = 0;
//        VertexID u;
//        for (auto v: vertices[cur]) {
//            if (max < deg[cur][v]) {
//                max = deg[cur][v];
//                u = v;
//            }
//        }
//
//        ui min = UINT_MAX;
//        for (auto v: adj[cur][u]) {
//            if (deg[1 - cur][v] * max <= w)
//                continue;
//            else {
//                if (min > deg[1 - cur][v] * max)
//                    min = deg[1 - cur][v] * max;
//            }
//        }

        w = min;
//        printf("w value: %d\n", w);
    }
//    printf("S: %d, T: %d\n", vertices[0].size(), vertices[1].size());
//    printf("w* value: %d\n", w);
    Graph w_core(true, graph.getVerticesCount());
    std::vector<bool> is_selected(graph.getVerticesCount(), false);

    for (auto u: vertices[1])
        is_selected[u] = true;

//#pragma omp parallel for
    for (auto u: vertices[0]) {
        for (auto v: graph.getOutNeighbors(u))
            if (is_selected[v])
                w_core.addDirectedEdge(u, v);
    }
    subgraph = w_core;
//    printf("%d\n", w_core.getEdgesCount());
}

void WCore::wCoreDecomposition(Graph &graph) {
    auto n = graph.getVerticesCount();
    auto adj = graph.getAdjList();
    std::vector<std::vector<ui>> deg(2);
//    std::vector<VertexID> vertices[2];
    deg[0] = graph.getOutDegrees();
    deg[1] = graph.getInDegrees();
    w = 1;
//    w = 27;
    for (ui i = 0; i < 2; i++)
        for (VertexID u = 0; u < n; u++)
            if (deg[i][u])
                vertices[i].push_back(u);
    printf("S: %d, T: %d\n", vertices[0].size(), vertices[1].size());

    bool changed;
    while (vertices[0].size() && vertices[1].size()) {
        changed = true;
//        std::vector<std::pair<VertexID, VertexID>> peeled_edges;
        while (changed) {
            changed = false;
            for (auto u: vertices[0]) {
                ui d = deg[0][u];
                ui dout = 0;
//                auto out = graph.getOutNeighbors(u);
                std::vector<VertexID> new_neighbors;
                for (auto v: adj[0][u]) {
                    if (d * deg[1][v] > w) {
                        dout++;
                        new_neighbors.push_back(v);
                    } else {
//                        peeled_edges.emplace_back(u, v);
                    }
                }
                if (dout < deg[0][u]) {
                    deg[0][u] = dout;
                    adj[0][u] = new_neighbors;
                    changed = true;
                }
            }

            for (auto v: vertices[1]) {
                ui d = deg[1][v];
                ui din = 0;
                std::vector<VertexID> new_neighbors;
                for (auto u: adj[1][v])
                    if (d * deg[0][u] > w) {
                        din++;
                        new_neighbors.push_back(u);
                    }
                if (din < deg[1][v]) {
                    deg[1][v] = din;
                    adj[1][v] = new_neighbors;
                    changed = true;
                }
            }
        }
//        for (auto edge: peeled_edges) {
//            edges.emplace_back(edge.first, edge.second);
//            induce_numbers.emplace_back(w);
//        }
        std::vector<ui> cnt(2, 0);
        for (ui i = 0; i < 2; i++)
            for (VertexID u = 0; u < n; u++)
                if (deg[i][u])
                    cnt[i]++;

        if (!cnt[0] || !cnt[1])
            break;

//        for (VertexID u = 0; u < n; u++)
//            if (deg[i][])

//        for (ui i = 0; i < 2; i++)
//            for(VertexID u = 0; u < n; u++)
//                if (!deg[i][u] && degrees[i][u])
//                    w_max.push_back(std::make_pair(w, ))

        degrees[0] = deg[0];
        degrees[1] = deg[1];
        vertices[0].clear();
        vertices[1].clear();

        for (ui i = 0; i < 2; i++)
            for (VertexID u = 0; u < n; u++)
                if (deg[i][u])
                    vertices[i].push_back(u);
        printf("S: %d, T: %d\n", vertices[0].size(), vertices[1].size());

        ui min = UINT_MAX;
        for (auto u: vertices[0]) {
            ui d = deg[0][u];
            for (auto v: adj[0][u])
                if (min > d * deg[1][v])
                    min = d * deg[1][v];
        }

//        ui cur = vertices[0].size() > vertices[1].size() ? 1 : 0;
//        ui max = 0;
//        VertexID u;
//        for (auto v: vertices[cur]) {
//            if (max < deg[cur][v]) {
//                max = deg[cur][v];
//                u = v;
//            }
//        }
//
//        ui min = UINT_MAX;
//        for (auto v: adj[cur][u]) {
//            if (deg[1 - cur][v] * max <= w)
//                continue;
//            else {
//                if (min > deg[1 - cur][v] * max)
//                    min = deg[1 - cur][v] * max;
//            }
//        }

        w = min;
        printf("w value: %d\n", w);
    }
    printf("S: %d, T: %d\n", vertices[0].size(), vertices[1].size());
    printf("w* value: %d\n", w);
}

void WCore::getMaxCNPair(Graph &graph, std::pair<ui, ui> &max_core_num_pair) {
//    omp_set_num_threads(32);
    auto n = graph.getVerticesCount();
    std::vector<ui> p(0);
//#pragma omp parallel for
    for (auto u: vertices[0]) {
        ui d = degrees[0][u];
        for (auto v: graph.getOutNeighbors(u))
            if (d * degrees[1][v] == w)
                p.push_back(degrees[1][v]);
    }
    sort(p.begin(), p.end());
    bool changed;
    ui d_cur = 0;
    ui pos = 0;
    max_core_num_pair.first = 0;
    max_core_num_pair.second = 0;
    while (!vertices[0].empty() && !vertices[1].empty()) {
        changed = true;
        while (d_cur == p[pos])
            pos++;
        if (pos >= p.size())
            break;
        d_cur = p[pos];
        while (changed) {
            changed = false;
            for (auto u: vertices[0]) {
                ui d = degrees[0][u];
                ui dout = 0;
                auto out = graph.getOutNeighbors(u);
                for (auto v: graph.getOutNeighbors(u)) {
                    if (d * degrees[1][v] > w || (degrees[1][v] != d_cur && d * degrees[1][v] == w))
                        dout++;
                    else if (degrees[1][v] == d_cur && d * degrees[1][v] == w)
                        max_core_num_pair.second = d_cur;
                }
                if (dout < degrees[0][u]) {
                    degrees[0][u] = dout;
                    changed = true;
                }
            }
            for (auto v: vertices[1]) {
                ui d = degrees[1][v];
                ui din = 0;
                for (auto u: graph.getInNeighbors(v)) {
                    if (d * degrees[0][u] > w || (d != d_cur && d * degrees[0][u] == w))
                        din++;
                    else if (d == d_cur && d_cur && d * degrees[0][u] == w)
                        max_core_num_pair.second = d_cur;
                }
                if (din < degrees[1][v]) {
                    degrees[1][v] = din;
                    changed = true;
                }
            }
        }

        for (ui i = 0; i < 2; i++)
            for (VertexID u = 0; u < n; u++)
                if (degrees[i][u])
                    vertices[i].push_back(u);
//        printf("S: %d, T: %d\n", vertices[0].size(), vertices[1].size());
    }
    max_core_num_pair.first = w / max_core_num_pair.second;
//    printf("x: %d, y: %d\n", max_core_num_pair.first, max_core_num_pair.second);
}

