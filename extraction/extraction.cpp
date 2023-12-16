//
// Created by yy on 12/1/23.
//

#include "extraction.h"

void Extraction::flowExactExtraction(Graph &graph, FlowNetwork &flow, double &l, double &r, std::vector<VertexID> *vertices) {
    std::vector<VertexID> tmp_S;
    ui n = graph.getVerticesCount();
    VertexID s = 0, t = 2 * n + 1;
    vertices[0].clear();
    vertices[1].clear();
    double mid = (l + r) / 2.0;

    flow.getMinCut(s, t, tmp_S);
    ui edge_num = 0;
    for (auto &v: tmp_S) {
        if (v == 0) continue;
        if (v <= n) {
            vertices[0].push_back(v - 1);
            for (auto &edge: flow.adj_[v]) {
                if (edge.to > n && edge.to <= 2 * n && flow.dist_[edge.to] >= n) {
                    edge_num++;
                }
            }
        } else vertices[1].push_back(v - n - 1);
    }

    if (!vertices[0].empty() && !vertices[1].empty()) {
        l = mid;
        if (graph.subgraph_density < edge_num / sqrt(vertices[0].size() * vertices[1].size())) {
            graph.subgraph_density = edge_num / sqrt(vertices[0].size() * vertices[1].size());
            graph.vertices[0] = vertices[0];
            graph.vertices[1] = vertices[1];
        }
    } else {
        r = mid;
    }
}

void Extraction::UndirectedflowExactExtraction(Graph &graph, FlowNetwork &flow, double &l, double &r, std::vector<VertexID> *vertices) {
    std::vector<VertexID> tmp_S;
    ui n = graph.getVerticesCount();
    VertexID s = 0, t = n + 1;
    vertices[0].clear();
    double mid = (l + r) / 2.0;
    flow.getMinCut(s, t, tmp_S);
    ui edge_num = 0;
    for (auto &v: tmp_S) {
        if (v == 0) continue;
        if (v <= n) {
            vertices[0].push_back(v - 1);
            for (auto &edge: flow.adj_[v]) {
                if (edge.to > 0 && edge.to <= n && flow.dist_[edge.to] >= n) {
                    edge_num++;
                }
            }
        }
    }
    if (!vertices[0].empty()) {
        l = mid;
        if (graph.subgraph_density < edge_num / vertices[0].size()) {
            graph.subgraph_density = edge_num / vertices[0].size();
            graph.vertices[0] = vertices[0];
        }
    } else {
        r = mid;
    }
}

void Extraction::UndirectedlpExactExtraction(Graph &graph, LinearProgamming &lp, std::vector<VertexID> *vertices) {
    std::vector<ui> y;
    std::vector<std::pair<double,VertexID>> tmp;
    ui n = graph.getVerticesCount();
    ui m = graph.getEdgesCount();
    y.resize(n);
    tmp.resize(n);
    vertices[0].clear();
    for(ui u = 0; u < n; u++){
        y[u] = 0;
        tmp[u] = std::make_pair(-lp.r[0][u],u);
    }
    for(ui i = 0; i < m; i++){
        if(lp.r[0][lp.alpha[0][i].id_first] < lp.r[0][lp.alpha[0][i].id_second] || 
            (lp.r[0][lp.alpha[0][i].id_first] == lp.r[0][lp.alpha[0][i].id_second] && lp.alpha[0][i].id_first > lp.alpha[0][i].id_second))
            y[lp.alpha[0][i].id_first]++;
        else
            y[lp.alpha[0][i].id_second]++;
    }
    sort(tmp.begin(), tmp.end());
    ui last_pos = 0;
    double sum = 0;
    double last_ans = 0;
    for(ui u = 0; u < n; u++){
        sum += y[tmp[u].second];
        if(sum / (u + 1) > last_ans){
            last_ans = sum / (u + 1);
            last_pos = u;
        }
    }
    for(ui u = 0; u <= last_pos; u++)
        vertices[0].push_back(tmp[u].second);
    graph.vertices[0] = vertices[0];
}