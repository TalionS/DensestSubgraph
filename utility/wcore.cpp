#include "wcore.h"

void WCore::wCoreDecomposition(Graph &graph) {
    auto n = graph.getVerticesCount();
    std::vector<std::vector<ui>> deg(2);
//    std::vector<VertexID> vertices[2];
    deg[0] = graph.getOutDegrees();
    deg[1] = graph.getInDegrees();
//    w = std::max(*std::max_element(deg[0].begin(), deg[0].end()),
//                     *std::max_element(deg[1].begin(), deg[1].end())) - 1;
    w = 27;
    for (ui i = 0; i < 2; i++)
        for (VertexID u = 0; u < n; u++)
            if (deg[i][u])
                vertices[i].push_back(u);
    printf("S: %d, T: %d\n", vertices[0].size(), vertices[1].size());

    bool changed;
    while (vertices[0].size() && vertices[1].size()) {
        changed = true;
        while (changed) {
            changed = false;
            for (auto u: vertices[0]) {
                ui d = deg[0][u];
                ui dout = 0;
                auto out = graph.getOutNeighbors(u);
                for (auto v: graph.getOutNeighbors(u))
                    if (d * deg[1][v] > w)
                        dout++;
                if (dout < deg[0][u]) {
                    deg[0][u] = dout;
                    changed = true;
                }
            }

            for (auto v: vertices[1]) {
                ui d = deg[1][v];
                ui din = 0;
                for (auto u: graph.getInNeighbors(v))
                    if (d * deg[0][u] > w)
                        din++;
                if (din < deg[1][v]) {
                    deg[1][v] = din;
                    changed = true;
                }
            }
        }
        std::vector<ui> cnt(2, 0);
        for (ui i = 0; i < 2; i++)
            for (VertexID u = 0; u < n; u++)
                if (deg[i][u])
                    cnt[i]++;

        if (!cnt[0] || !cnt[1])
            break;

        degrees[0] = deg[0];
        degrees[1] = deg[1];
        vertices[0].clear();
        vertices[1].clear();

        for (ui i = 0; i < 2; i++)
            for (VertexID u = 0; u < n; u++)
                if (deg[i][u])
                    vertices[i].push_back(u);
        printf("S: %d, T: %d\n", vertices[0].size(), vertices[1].size());

        ui cur = vertices[0].size() > vertices[1].size() ? 1 : 0;
        ui max = 0;
        VertexID u;
        for (auto v: vertices[cur]) {
            if (max < deg[cur][v]) {
                max = deg[cur][v];
                u = v;
            }
        }

        ui min = UINT_MAX;
        for (auto v: cur ? graph.getInNeighbors(u) : graph.getOutNeighbors(u)) {
            if (deg[1 - cur][v] * max <= w)
                continue;
            else {
                if (min > deg[1 - cur][v] * max)
                    min = deg[1 - cur][v] * max;
            }
        }

        w = min;
        printf("w value: %d\n", w);
    }
    printf("S: %d, T: %d\n", vertices[0].size(), vertices[1].size());
    printf("w* value: %d\n", (int) w);
}

void WCore::getMaxCNPair(Graph &graph, std::pair<ui, ui> &max_core_num_pair) {
    auto n = graph.getVerticesCount();
    std::vector<ui> p(0);
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
    while (vertices[0].size() && vertices[1].size()) {
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
        printf("S: %d, T: %d\n", vertices[0].size(), vertices[1].size());
    }
    max_core_num_pair.first = w / max_core_num_pair.second;
    printf("x: %d, y: %d\n", max_core_num_pair.first, max_core_num_pair.second);
}

