//
// Created by yy on 11/23/23.
//
#include "reduction.h"
#include "utility/graph.h"
#include <algorithm>

std::vector<ui> Reduction::coreDecomposition(const Graph &graph) {
    ui vertices_count = graph.getVerticesCount();
    std::vector<ui> degrees = graph.getDegrees();
    std::vector<ui> core(vertices_count, 0);
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
    return core;
};

std::vector<VertexID> *Reduction::xyCoreDecomposition(const Graph &graph, ui x, ui y) {
    auto vertices_count = graph.getVerticesCount();
    auto *degrees = new std::vector<ui>[2];
    degrees[0] = graph.getOutDegrees();
    degrees[1] = graph.getInDegrees();
    auto *vert = new std::vector<VertexID>[2];
    auto *bin = new std::vector<ui>[2];
    auto *pos = new std::vector<ui>[2];
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

    VertexID i = 0, j = 0;
    while(degrees[0][vert[0][i]] < x || degrees[1][vert[1][j]] < y) {
        for (; degrees[0][vert[0][i]] < x; i++) {
            auto u = vert[0][i];
            for (auto v : graph.getOutNeighbors(u)) {
                if (degrees[1][v] >= y) {
//                    dec_deg(1, v);
                    auto dt = degrees[1][v], pt = pos[1][v];
                    auto pw = bin[1][dt], w = vert[1][pw];
                    if (v != w) {
                        pos[1][v] = pw; vert[1][pt] = w;
                        pos[1][w] = pt; vert[1][pw] = v;
                    }
                    ++bin[1][dt];
                    --degrees[1][v];
                }
            }
        }

        for (; degrees[1][vert[1][j]] < y; j++) {
            auto u = vert[1][j];
            for (auto v : graph.getInNeighbors(u)) {
                if (degrees[0][v] >= x) {
//                    dec_deg(0, v);
                    auto dt = degrees[0][v], pt = pos[0][v];
                    auto pw = bin[0][dt], w = vert[0][pw];
                    if (v != w) {
                        pos[0][v] = pw; vert[0][pt] = w;
                        pos[0][w] = pt; vert[0][pw] = v;
                    }
                    ++bin[0][dt];
                    --degrees[0][v];
                }
            }
        }
    }
    auto *vertices = new std::vector<VertexID>[2];
    for (; i < vertices_count; i++) {
        auto u = vert[0][i];
        vertices[0].push_back(u);
        for (auto v : graph.getOutNeighbors(u)) {
            if (degrees[1][v] >= y) {
                vertices[1].push_back(v);
            }
        }
    }
    sort( vertices[1].begin(), vertices[1].end() );
    vertices[1].erase( unique( vertices[1].begin(), vertices[1].end() ), vertices[1].end() );
    return vertices;
}