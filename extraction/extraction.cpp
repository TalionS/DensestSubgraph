//
// Created by yy on 12/1/23.
//

#include "extraction.h"
#include "wcore.h"

void
Extraction::flowExactExtraction(Graph &graph, Graph &subgraph, std::pair<double, double> ratio, FlowNetwork &flow,
                                double &l, double &r, double &ratio_o, double &ratio_p, bool is_map) {
    std::vector<VertexID> tmp_S;
    ui n = subgraph.getVerticesCount();
    std::vector<std::vector<VertexID>> vertices(2);
    double mid = (l + r) / 2.0;
//    std::vector<ui> degrees[2];
//    degrees[1] = graph.getInDegrees();
//    degrees[0] = graph.getOutDegrees();
//    std::vector<ui> map(1, 0);
//    std::vector<ui> cnt(2, 0);
//    ui tmp = 1;
//    for (ui i = 0; i < 2; i++){
//        for (VertexID u = 0; u < n; u++){
//            if (degrees[i][u]) {
//                map.push_back(u);
//                tmp++;
//            }
//        }
//        cnt[i] = tmp;
//    }
//    VertexID s = 0, t = tmp;
//    flow.getMinCut(s, t, tmp_S);
//    ui edge_num = 0;
//    for (auto &v: tmp_S) {
//        if (v == 0) continue;
//        if (v < cnt[0]) {
//            vertices[0].push_back(map[v]);
//            for (auto &edge: flow.adj_[v]) {
//                if (edge.to >= cnt[0] && edge.to < cnt[1] && flow.dist_[edge.to] >= tmp + 1) {
//                    edge_num++;
//                }
//            }
//        } else vertices[1].push_back(map[v]);
//    }
    ui s = 0, t = 2 * n + 1;
    flow.getMinCut(s, t, tmp_S);
    ui edge_num = 0;
    if (!subgraph.is_mapped){
        for (auto &v: tmp_S) {
            if (v == 0) continue;
            if (v <= n) {
                vertices[0].push_back(v - 1);
                for (auto &edge: flow.adj_[v]) {
                    if (edge.to > n && edge.to <= 2 * n && flow.dist_[edge.to] > t) {
                        edge_num++;
                    }
                }
            } else vertices[1].push_back(v - n - 1);
        }
    } else {
        for (auto &v: tmp_S) {
            if (v == 0) continue;
            if (v <= n) {
                vertices[0].push_back(subgraph.map[v - 1]);
                for (auto &edge: flow.adj_[v]) {
                    if (edge.to > n && edge.to <= 2 * n && flow.dist_[edge.to] > t) {
                        edge_num++;
                    }
                }
            } else vertices[1].push_back(subgraph.map[v - n - 1]);
        }
    }

    if (!vertices[0].empty() && !vertices[1].empty()) {
        l = mid;
        double ratios;
        if (is_map) {
            if (ratio.first < 1 && ratio.second > 1) {
                ratios = 1;
            } else if (ratio.second <= 1) {
                ratios = (ratio.first + ratio.second) / 2;
            } else if (ratio.first >= 1) {
                ratios = 2 / (1 / ratio.first + 1 / ratio.second);
            }
        } else
            ratios = (ratio.first + ratio.second) / 2;
        ratio_o = (double) vertices[0].size() / vertices[1].size();
        ratio_p = ratios * ratios / ratio_o;
        if (ratio_o > ratio_p)
            std::swap(ratio_o, ratio_p);
        if (graph.subgraph_density < edge_num / sqrt((double) vertices[0].size() * vertices[1].size())) {
            graph.subgraph_density = edge_num / sqrt((double) vertices[0].size() * vertices[1].size());
            graph.vertices[0] = vertices[0];
            graph.vertices[1] = vertices[1];
        }
    } else {
        r = mid;
    }
}

void Extraction::directedCoreApproExtraction(Graph &graph, Graph &subgraph, std::pair<ui, ui> &max_core_num_pair) {
//    return;
    XYCore xy_core;
    xy_core.xyCoreInitialization(subgraph, true);
    Graph max_xy_core(true, 0);
    xy_core.generateXYCore(subgraph, max_xy_core, max_core_num_pair.first, max_core_num_pair.second, false, false,
                           false);
    std::vector<ui> degrees[2];
    degrees[0] = max_xy_core.getOutDegrees();
    degrees[1] = max_xy_core.getInDegrees();
    std::vector<VertexID> vertices[2];
    for (ui i = 0; i < 2; i++) {
        for (VertexID u = 0; u < graph.getVerticesCount(); u++) {
            if (degrees[i][u])
                vertices[i].emplace_back(u);
        }
    }
    graph.subgraph_density = max_xy_core.getEdgesCount() / sqrt((double) vertices[0].size() * vertices[1].size());
    graph.vertices[0] = vertices[0];
    graph.vertices[1] = vertices[1];
}

void Extraction::directedBSApproExtraction(Graph &graph, std::vector<std::vector<bool>> is_peeled,
                                           std::vector<std::vector<VertexID>> &vertices) {
    for (ui i = 0; i < 2; i++) {
        vertices[i].clear();
        for (ui u = 0; u < graph.getVerticesCount(); u++)
            if (!is_peeled[i][u])
                vertices[i].push_back(u);
    }

}

void Extraction::directedPMApproExtraction(Graph &graph, ui edges_count, std::vector<std::vector<VertexID>> vertices) {
    if (!vertices[0].size() || !vertices[1].size())
        return;
    if (edges_count > graph.subgraph_density * sqrt(vertices[0].size() * vertices[1].size())) {
        graph.subgraph_density = edges_count / sqrt(vertices[0].size() * vertices[1].size());
        graph.vertices[0] = vertices[0];
        graph.vertices[1] = vertices[1];
    }
}

void Extraction::directedCPExtraction(Graph &graph, LinearProgramming &lp, std::pair<ui, ui> &best_pos,
                                      std::vector<std::vector<VertexID>> &vertices, std::pair<double, double> ratios,
                                      double &ratio_o, double &ratio_p, double &rho, double &rho_c, bool is_map) {
    if (!graph.getEdgesCount())
        return;
    rho = 0;
    rho_c = 0;
    best_pos = std::make_pair(0, 0);
    double ratio;
    if (is_map) {
        if (ratios.first < 1 && ratios.second > 1) {
            ratio = 1;
        } else if (ratios.second <= 1) {
            ratio = (ratios.first + ratios.second) / 2;
        } else if (ratios.first >= 1) {
            ratio = 2 / (1 / ratios.first + 1 / ratios.second);
        }
    } else
        ratio = (ratios.first + ratios.second) / 2;
    ratio_o = 0;
    ratio_p = 0;
//    double ratio = ratios.first / ratios.second;
    auto out_degrees = graph.getOutDegrees();
    auto in_degrees = graph.getInDegrees();
    std::vector<std::vector<ui>> y(2);
    std::vector<std::vector<std::pair<double, VertexID>>> tmp_r(2);
    std::vector<ui> cnt(2, 0);
    ui n = graph.getVerticesCount();
//    for (ui i = 0; i < 2; i++)
//        y[i].resize(n);
    for (VertexID u = 0; u < n; u++) {
        if (out_degrees[u]) {
            cnt[0]++;
            tmp_r[0].emplace_back(std::make_pair(-lp.r[0][u], u));
        }
        if (in_degrees[u]) {
            cnt[1]++;
            tmp_r[1].emplace_back(std::make_pair(-lp.r[1][u], u));
        }
    }
//    cnt[0] = graph.vertex_ids[0].size();
//    cnt[1] = graph.vertex_ids[1].size();
//    for (ui i = 0; i < 2; i++)
//        for (auto &u: graph.vertex_ids[i])
//            tmp_r[i].emplace_back(std::make_pair(-lp.r[i][u], u));

    for (ui i = 0; i < 2; i++)
        sort(tmp_r[i].begin(), tmp_r[i].end());
    double sum = 0;
    std::vector<double> pos(2, 0);
    best_pos.first = tmp_r[0][0] < tmp_r[1][0] ? 0 : 1;
//    bool flag = false;
    std::vector<std::vector<bool>> is_selected(2);
    for (ui i = 0; i < 2; i++)
        is_selected[i].resize(n, false);
    while (pos[0] < cnt[0] || pos[1] < cnt[1]) {
        ui cur;
        if (pos[1] == cnt[1] || (pos[0] != cnt[0] && tmp_r[0][pos[0]] < tmp_r[1][pos[1]])) {
            cur = 0;
        } else {
            cur = 1;
        }
//        printf("%f\n", tmp_r[cur][pos_[cur]].first);
//        sum += y[cur][tmp_r[cur][pos_[cur]].second];
        is_selected[cur][tmp_r[cur][pos[cur]].second] = true;
        for (auto v: cur ? graph.getInNeighbors(tmp_r[cur][pos[cur]].second) : graph.getOutNeighbors(
                tmp_r[cur][pos[cur]].second)) {
            if (is_selected[1 - cur][v])
                sum++;
        }
        pos[cur]++;
        if (pos[0] == 0 || pos[1] == 0) continue;
        double ratio_prime = pos[0] / pos[1];
//        if (!epsilon){
        if (2 * sqrt(ratio * ratio_prime) * sum > rho_c * (ratio + ratio_prime) * sqrt(pos[0] * pos[1])) {
            best_pos = std::make_pair(cur, pos[cur]);
            ratio_o = ratio_prime;
            rho = sum / sqrt(pos[0] * pos[1]);
            rho_c = 2 * sqrt(ratio * ratio_prime) / (ratio + ratio_prime) * sum / sqrt(pos[0] * pos[1]);
        }
    }
//        } else {
//            if (sum / sqrt(pos[0] * pos[1]) > rho) {
//                rho = sum / sqrt(pos[0] * pos[1]);
//                best_pos = std::make_pair(cur, pos[cur]);
//            }
//            if (2 * sqrt(ratio * ratio_prime) * sum > rho_c * (ratio + ratio_prime) * sqrt(pos[0] * pos[1])) {
//                ratio_o = ratio_prime;
//                rho_c = 2 * sqrt(ratio * ratio_prime) / (ratio + ratio_prime) * sum / sqrt(pos[0] * pos[1]);
//            }
//        }
//    }
    ratio_p = ratio * ratio / ratio_o;
    if (ratio_o > ratio_p) std::swap(ratio_o, ratio_p);
    vertices[0].clear();
    vertices[1].clear();
    pos.assign(2, 0);
    ui cur = tmp_r[0][0] < tmp_r[1][0] ? 0 : 1;
    while (cur != best_pos.first || pos[cur] != best_pos.second) {
        if (pos[1] == cnt[1] || (pos[0] != cnt[0] && tmp_r[0][pos[0]] < tmp_r[1][pos[1]])) {
            cur = 0;
        } else {
            cur = 1;
        }
        vertices[cur].emplace_back(tmp_r[cur][pos[cur]].second);
        pos[cur]++;
    }
}

//void Extraction::directedVWApproExtraction(Graph &graph, LinearProgramming &lp,
//                                           std::vector<std::vector<VertexID>> &vertices,
//                                           std::pair<double, double> ratios, double &rho, double &vw_rho) {
//    rho = 0;
//    vw_rho = 0;
//    std::pair<ui, ui> best_pos = std::make_pair(0, 0);
//    double ratio = ratios.first / ratios.second;
////    printf("%f\n", sqrt(1 / graph.getVerticesCount()));
////    double ratio = ratios.first / ratios.second;
//    auto out_degrees = graph.getOutDegrees();
//    auto in_degrees = graph.getInDegrees();
//    std::vector<std::vector<ui>> y(2);
//    std::vector<std::vector<std::pair<double, VertexID>>> tmp_r(2);
//    std::vector<ui> cnt(2, 0);
//    ui n = graph.getVerticesCount();
//    for (ui i = 0; i < 2; i++)
//        y[i].resize(n);
//    for(VertexID u = 0; u < n; u++){
//        if (out_degrees[u]){
//            cnt[0]++;
//            tmp_r[0].emplace_back(std::make_pair(-lp.r[0][u], u));
//        }
//        if (in_degrees[u]){
//            cnt[1]++;
//            tmp_r[1].emplace_back(std::make_pair(-lp.r[1][u], u));
//        }
//    }
//    for (ui i = 0; i < 2; i++)
//        sort(tmp_r[i].begin(), tmp_r[i].end());
//    double sum = 0;
//    std::vector<ui> pos(2, 0);
//    best_pos.first = tmp_r[0][0] > tmp_r[1][0] ? 1 : 0;
//    std::vector<std::vector<bool>> is_selected(2);
//    for (ui i = 0; i < 2; i++)
//        is_selected[i].resize(n, false);
//    while (pos[0] < cnt[0] || pos[1] < cnt[1]){
//        ui cur;
//        if (pos[1] == cnt[1] || (pos[0] != cnt[0] && tmp_r[0][pos[0]] < tmp_r[1][pos[1]])){
//            cur = 0;
//        }
//        else{
//            cur = 1;
//        }
////        printf("%f\n", tmp_r[cur][pos_[cur]].first);
////        sum += y[cur][tmp_r[cur][pos_[cur]].second];
//        is_selected[cur][tmp_r[cur][pos[cur]].second] = true;
//        for (auto v: cur? graph.getInNeighbors(tmp_r[cur][pos[cur]].second) : graph.getOutNeighbors(tmp_r[cur][pos[cur]].second)){
//            if (is_selected[1 - cur][v])
//                sum++;
//        }
//        pos[cur]++;
//        if (pos[0] == 0 || pos[1] == 0) continue;
////        vw_rho = 2 * edge_num / (cnt[0] / sqrt_c + cnt[1] * sqrt_c)
////        printf("vw_rho: %.20f\n", (double) 2 * sum / (pos_[0] / sqrt(ratio) + pos_[1] * sqrt(ratio)));
//        if ((double) 2 * sum / (pos[0] / sqrt(ratio) + pos[1] * sqrt(ratio)) > vw_rho){
//            best_pos = std::make_pair(cur, pos[cur]);
//            rho = sum / sqrt(pos[0] * pos[1]);
//            vw_rho = (double) 2 * sum / (pos[0] / sqrt(ratio) + pos[1] * sqrt(ratio));
//        }
//    }
//    vertices[0].clear();
//    vertices[1].clear();
//    pos.assign(2, 0);
//    ui cur = tmp_r[0][0] > tmp_r[1][0] ? 1 : 0;
//    while (cur != best_pos.first || pos[cur] != best_pos.second){
//        if (pos[1] == cnt[1] || (pos[0] != cnt[0] && tmp_r[0][pos[0]] < tmp_r[1][pos[1]])){
//            cur = 0;
//        }
//        else{
//            cur = 1;
//        }
//        vertices[cur].emplace_back(tmp_r[cur][pos[cur]].second);
//        pos[cur]++;
//    }
//}

void Extraction::directedVWApproExtraction(Graph &graph, Graph &vw_graph, std::vector<std::vector<VertexID>> vertices,
                                           std::vector<VertexID> *verticess, double &ratio_o, double &ratio_p) {
    auto n = graph.getVerticesCount();
    vertices[0].clear();
    vertices[1].clear();
    std::vector<bool> is_selected(n, false);
    for (auto u: verticess[0]) {
        if (vw_graph.map[u] < n) {
            vertices[0].emplace_back(u);
            is_selected[u] = true;
        } else
            vertices[1].emplace_back(u);
    }
    ratio_o = (double) vertices[0].size() / vertices[1].size();
    n = vw_graph.getVerticesCount() - 1;
    ratio_p = vw_graph.weight_[n] * vw_graph.weight_[n] * vw_graph.weight_[n] * vw_graph.weight_[n] * 16 / ratio_o;
    if (ratio_o > ratio_p)
        std::swap(ratio_o, ratio_p);
    //std::cout<<vertices[0].size()<<" "<<vertices[1].size()<<std::endl;
    double edges_count = 0;
    for (auto v: vertices[1])
        for (auto u: vw_graph.getNeighbors(v))
            if (is_selected[u])
                edges_count++;
    if (edges_count / sqrt((double) vertices[0].size() * vertices[1].size()) > graph.subgraph_density) {
        graph.subgraph_density = edges_count / sqrt((double) vertices[0].size() * vertices[1].size());
        graph.vertices[0] = vertices[0];
        graph.vertices[1] = vertices[1];
    }
    //std::cout<<"haha"<<std::endl;
}

void Extraction::UndirectedflowExactExtraction(Graph &graph, FlowNetwork &flow, double &l, double &r,
                                               std::vector<VertexID> *vertices) {
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
    std::vector<std::pair<double, VertexID>> tmp;
    ui n = graph.getVerticesCount();
    ui m = graph.getEdgesCount();
    y.resize(n);
    tmp.resize(n);
    vertices[0].clear();
    for (ui u = 0; u < n; u++) {
        y[u] = 0;
        tmp[u] = std::make_pair(-lp.r[0][u], u);
    }
    for (ui i = 0; i < m; i++) {
        if (lp.r[0][lp.alpha[i].id_first] < lp.r[0][lp.alpha[i].id_second] ||
            (lp.r[0][lp.alpha[i].id_first] == lp.r[0][lp.alpha[i].id_second] &&
             lp.alpha[i].id_first > lp.alpha[i].id_second))
            y[lp.alpha[i].id_first] += 1.0;
        else
            y[lp.alpha[i].id_second] += 1.0;
    }
    sort(tmp.begin(), tmp.end());
    ui last_pos = 0;
    double sum = 0, opt = 0, weight_sum = 0;
    for(ui i = 0; i < n; i++){
        sum += y[tmp[i].second];
        weight_sum += lp.weight[tmp[i].second];
        if(sum / weight_sum > opt){
            opt = sum / weight_sum;
            last_pos = i + 1;
        }
    }
    graph.subgraph_density_lower_bound = std::max(graph.subgraph_density_lower_bound, opt);
    for (ui u = 0; u < last_pos; u++)
        vertices[0].push_back(tmp[u].second);
    graph.vertices[0] = vertices[0];
}