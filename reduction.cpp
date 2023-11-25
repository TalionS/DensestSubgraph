//
// Created by yy on 11/23/23.
//
#include "reduction.h"
#include "graph.h"
#include <algorithm>

std::vector<ui> Reduction::coreDecomposition(const Graph& graph){
    ui vertices_count = graph.getVerticesCount();
    std::vector<ui> degrees = graph.getDegrees();
    std::vector<ui> core(vertices_count, 0);
    std::vector<ui> rid(vertices_count, 0);
    std::vector<ui> id(vertices_count, 0);
    for(ui i = 0; i < vertices_count; i++) ++id[degrees[i]];
    for(ui i = 1; i < vertices_count; i++) id[i] += id[i-1];

    for(ui i = 0; i < vertices_count; i++) rid[i] = --id[degrees[i]];
    for(ui i = 0; i < vertices_count; i++) id[rid[i]] = i;

    std::vector<ui> degree_start(vertices_count+1);
    for(ui i = 0, j = 0; i <= vertices_count; i++){
        while(j < vertices_count && degrees[id[j]] < i) ++j;
        degree_start[i] = j;
    }

    ui max_core = 0;
    for(ui i = 0; i < vertices_count; i++){
        ui u = id[i];
//        assert(degree_start[degrees[u]] == i);
        if(degrees[u] > max_core) max_core = degrees[u];
        core[u] = max_core;

        ++degree_start[degrees[u]];
        if(degrees[u] == 0) continue;

        degree_start[degrees[u] - 1] = degree_start[degrees[u]];
        for(ui j: graph.getNeighbors(u)) if(rid[j] > i){
            ui pos1 = degree_start[degrees[j]], pos2 = rid[j];
            std::swap(id[pos1], id[pos2]);
            rid[id[pos1]] = pos1; rid[id[pos2]] = pos2;
            ++degree_start[degrees[j]];
            --degrees[j];
        }
    }
    return core;
};