//
// Created by yy on 12/1/23.
//

#include "extraction.h"

void Extraction::flowExactExtraction(Graph &graph, FlowNetwork &flow, double &l, double &r, ui &s_size, ui &t_size) {
    std::vector<VertexID> tmp_S;
    ui n = graph.getVerticesCount();
    VertexID s = 0, t = 2 * n + 1;
    std::vector<std::vector<VertexID>> vertices(2);
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
        s_size = vertices[0].size();
        t_size = vertices[1].size();
        if (graph.subgraph_density < edge_num / sqrt(vertices[0].size() * vertices[1].size())) {
            graph.subgraph_density = edge_num / sqrt(vertices[0].size() * vertices[1].size());
            graph.vertices[0] = vertices[0];
            graph.vertices[1] = vertices[1];
        }
    } else {
        r = mid;
    }
}

void Extraction::directedLPExactExtraction(Graph &graph,
                                           LinearProgramming &lp,
                                           std::pair<ui, ui> &best_pos,
                                           std::vector<std::vector<VertexID>> &vertices,
                                           std::pair<double, double> ratios,
                                           double &ratio_o,
                                           double &ratio_p,
                                           double &rho_c) {
    rho_c = 0;
    ratio_o = 0;
    ratio_p = 0;
    best_pos = std::make_pair(0, 0);
    double ratio = (ratios.first + ratios.second) / 2;
//    double ratio = ratios.first / ratios.second;
    auto out_degrees = graph.getOutDegrees();
    auto in_degrees = graph.getInDegrees();
    std::vector<std::vector<ui>> y(2);
    std::vector<std::vector<std::pair<double, VertexID>>> tmp_r(2);
    std::vector<ui> cnt(2, 0);
    ui n = graph.getVerticesCount();
    ui m = graph.getEdgesCount();
    for (ui i = 0; i < 2; i++)
        y[i].resize(n);
    for(VertexID u = 0; u < n; u++){
        if (out_degrees[u]){
            cnt[0]++;
            tmp_r[0].emplace_back(std::make_pair(-lp.r[0][u], u));
        }
        if (in_degrees[u]){
            cnt[1]++;
            tmp_r[1].emplace_back(std::make_pair(-lp.r[1][u], u));
        }
    }
    for (ui i = 0; i < m; i++){
        if (lp.r[0][lp.alpha[0][i].id_first] < lp.r[1][lp.alpha[0][i].id_second] ||
                (lp.r[0][lp.alpha[0][i].id_first] == lp.r[1][lp.alpha[0][i].id_second] && lp.alpha[0][i].id_first > lp.alpha[0][i].id_second))
            y[0][lp.alpha[0][i].id_first]++;
        else
            y[1][lp.alpha[0][i].id_second]++;
    }
    for (ui i = 0; i < 2; i++)
        sort(tmp_r[i].begin(), tmp_r[i].end());
    double sum = 0;
    std::vector<ui> pos(2, 0);

    while (pos[0] < cnt[0] || pos[1] < cnt[1]){
        ui cur;
        if (pos[1] == cnt[1] || tmp_r[0][pos[0]] < tmp_r[1][pos[1]]){
            cur = 0;
        }
        else{
            cur = 1;
        }
//        printf("%f\n", tmp_r[cur][pos[cur]].first);
        sum += y[cur][tmp_r[cur][pos[cur]].second];
        pos[cur]++;
        if (pos[0] == 0 || pos[1] == 0) continue;
        double ratio_prime = pos[0] / pos[1];
        if (2 * sqrt(ratio * ratio_prime) / (ratio + ratio_prime) * sum / sqrt(pos[0] * pos[1]) > rho_c){
            best_pos = std::make_pair(cur, pos[cur]);
            ratio_o = ratio_prime;
            rho_c = 2 * sqrt(ratio * ratio_prime) / (ratio + ratio_prime) * sum / sqrt(pos[0] * pos[1]);
        }
    }
    if (ratio_o == 0) ratio_o = ratio;
    ratio_p = ratio * ratio / ratio_o;
    if (ratio_o > ratio_p) std::swap(ratio_o, ratio_p);
    vertices[0].clear();
    vertices[1].clear();
    pos.assign(2, 0);
    if (cnt[0] && cnt[1]){
        ui cur = tmp_r[0][0] > tmp_r[1][0] ? 1 : 0;
        while (cur != best_pos.first || pos[cur] != best_pos.second){
            if (pos[1] == cnt[1] || tmp_r[0][pos[0]] < tmp_r[1][pos[1]]){
                cur = 0;
            }
            else{
                cur = 1;
            }
            vertices[cur].emplace_back(tmp_r[cur][pos[cur]].second);
            pos[cur]++;
        }
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

void Extraction::UndirectedlpExactExtraction(Graph &graph, LinearProgramming &lp, std::vector<VertexID> *vertices) {
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