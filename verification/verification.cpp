//
// Created by yy on 12/1/23.
//

#include "verification.h"
#include <bits/stdc++.h>

bool Verification::flowExactVerification(Graph &graph, double l, double r) {
    ui n = graph.getVerticesCount();
    double bias = 1.0 / sqrt((double) n * (n - 1)) - 1.0 / n;

//    double bias = 1 / sqrt((double) (n - xy_core.bin_[0][1]) * (n - xy_core.bin_[1][1]) - std::min(n - xy_core.bin_[0][1], n - xy_core.bin_[1][1]))
//    - 1 / sqrt((double) (n - xy_core.bin_[0][1]) * (n - xy_core.bin_[1][1]));
    if (r - l > bias) return true;
    else return false;
}

bool Verification::directedBSApproVerification(Graph &graph, ui edges_count,
                                               std::vector<std::vector<VertexID>> vertices) {
    if (edges_count) {
        if (graph.subgraph_density * sqrt((double) vertices[0].size() * vertices[1].size()) < edges_count) {
            graph.subgraph_density = edges_count / sqrt((double) vertices[0].size() * vertices[1].size());
            graph.vertices[0] = vertices[0];
            graph.vertices[1] = vertices[1];
        }
    }
    return edges_count;
}

bool Verification::directedFixedKSApproVerification(Graph &graph, ui &cnt, ui edges_count,
                                                    std::vector<std::vector<VertexID>> vertices) {
    cnt++;
    if (edges_count) {
        if (graph.subgraph_density * sqrt((double) vertices[0].size() * vertices[1].size()) < edges_count) {
            graph.subgraph_density = edges_count / sqrt((double) vertices[0].size() * vertices[1].size());
            graph.vertices[0] = vertices[0];
            graph.vertices[1] = vertices[1];
        }
    }
    return cnt < 2;
}

bool Verification::directedPMApproVerification(std::vector<std::vector<VertexID>> vertices) {
    return vertices[0].size() && vertices[1].size();
}

bool Verification::UndirectedflowExactVerification(Graph &graph, double l, double r) {
    ui n = graph.getVerticesCount();
    double bias = 1.0 / ((double) n * (n - 1));
    if (r - l > bias) return true;
    else return false;
}

bool
Verification::directedCPVerification(Graph &graph, Graph &subgraph, LinearProgramming &lp, std::pair<ui, ui> best_pos,
                                     std::vector<std::vector<VertexID>> &vertices, std::pair<double, double> ratios,
                                     double rho, double rho_c, double &ratio_o, double &ratio_p,
                                     bool &stable_set_reduction, std::vector<std::pair<VertexID, VertexID>> &edges,
                                     double epsilon, bool stable, bool is_map) {
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
//    double ratio = ratios.first / ratios.second;
    stable_set_reduction = false;
    std::vector<std::vector<std::pair<double, VertexID>>> tmp_r(2);
    std::vector<ui> cnt(2, 0);
    auto out_degrees = subgraph.getOutDegrees();
    auto in_degrees = subgraph.getInDegrees();
    ui n = subgraph.getVerticesCount();
    ui m = subgraph.getEdgesCount();
    edges.clear();
    if (!m) {
        ratio_o = (ratios.first + ratios.second) / 2;
        ratio_p = ratio_o;
        return false;
    }

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
//    if (!cnt[0] || !cnt[1])
//        return true;
    for (ui i = 0; i < 2; i++)
        sort(tmp_r[i].begin(), tmp_r[i].end());

    ui cur = tmp_r[0][0] < tmp_r[1][0] ? 0 : 1;
    bool flag = true;
    if (-tmp_r[cur][0].first < graph.subgraph_density * sqrt(1 + epsilon)) {
        double t = -tmp_r[cur][0].first / graph.subgraph_density / sqrt(1 + epsilon);
        ratio_o = (2 * ratio - ratio * t * t - 2 * sqrt(ratio * ratio - ratio * ratio * t * t)) / (t * t);
        ratio_p = ratio * ratio / ratio_o;
        if (ratio_o > ratio_p)
            std::swap(ratio_o, ratio_p);
//        printf("early stop, ratio_o: %f, ratio_p: %f\n", ratio_o, ratio_p);
//        return false;
        flag = false;
    }
    if (flag){
        if (epsilon) {
//        flag = -tmp_r[cur][0].first / rho_c > sqrt(1 + epsilon);
//        flag = -tmp_r[cur][0].first > sqrt(1 + epsilon) * rho_c;
//        ratio_o = std::min(ratio_o, ratio / (1 + epsilon));
//        ratio_p = std::max(ratio_p, ratio * (1 + epsilon));
            if (-tmp_r[cur][0].first / rho_c <= sqrt(1 + epsilon)) {
                flag = false;
                ratio_o = std::min(ratio_o, ratio / (1 + epsilon));
                ratio_p = std::max(ratio_p, ratio * (1 + epsilon));
            } else if (-tmp_r[cur][0].first / rho_c <= 1 + epsilon && ratio_o < ratio / (1 + epsilon) &&
                       ratio_p > ratio * (1 + epsilon)) {
                flag = false;
            } else
                flag = true;
            if (stable) {
                std::vector<std::vector<bool>> selected(2);
                for (ui i = 0; i < 2; i++) {
                    selected[i].resize(n, false);
                    for (auto u: vertices[i])
                        selected[i][u] = true;
                }
                std::vector<double> r[2];
                r[0] = lp.r[0];
                r[1] = lp.r[1];
                std::vector<Alpha> alpha = lp.alpha;
                for (ui i = 0; i < lp.edges_count_; i++) {
//            if (!lp.alpha[i].is_selected)
//                continue;
                    if (selected[0][alpha[i].id_first] && !selected[1][alpha[i].id_second]) {
                        r[0][alpha[i].id_first] -= 2.0 * sqrt(ratio) * alpha[i].weight_first;
                        r[1][alpha[i].id_second] += 2.0 / sqrt(ratio) * alpha[i].weight_first;
                        alpha[i].weight_first = 0;
                        alpha[i].weight_second = 1;
//                    out_degrees[alpha[i].id_first]--;
                    } else if (!selected[0][alpha[i].id_first] && selected[1][alpha[i].id_second]) {
                        r[0][alpha[i].id_first] += 2.0 * sqrt(ratio) * alpha[i].weight_second;
                        r[1][alpha[i].id_second] -= 2.0 / sqrt(ratio) * alpha[i].weight_second;
                        alpha[i].weight_first = 1;
                        alpha[i].weight_second = 0;
                        in_degrees[alpha[i].id_second]--;
                    } else if (selected[0][alpha[i].id_first] && selected[1][alpha[i].id_second]) {
                        edges.emplace_back(alpha[i].id_first, alpha[i].id_second);
                    }
                }
                std::vector<ui> pos(2, 0);
                double min_r = r[cur][tmp_r[cur][0].second];
//        printf("(%f, %f)\n", min_r, rho_c);
                while (cur != best_pos.first || pos[cur] != best_pos.second) {
                    if (pos[1] == cnt[1] || (pos[0] != cnt[0] && tmp_r[0][pos[0]] < tmp_r[1][pos[1]])) {
                        cur = 0;
                    } else {
                        cur = 1;
                    }
//        min_r = std::max(min_r, tmp_r[cur][pos_[cur]].first);
                    min_r = std::min(min_r, r[cur][tmp_r[cur][pos[cur]].second]);
                    pos[cur]++;
                }
                while (pos[0] < cnt[0] && pos[1] < cnt[1]) {
                    if (pos[1] == cnt[1] || (pos[0] != cnt[0] && tmp_r[0][pos[0]] < tmp_r[1][pos[1]])) {
                        cur = 0;
                    } else {
                        cur = 1;
                    }
//        if (min_r > tmp_r[cur][pos_[cur]].first)
                    if (min_r < r[cur][tmp_r[cur][pos[cur]].second])
                        return true;
                    pos[cur]++;
                }
                if (edges.size() < m) {
                    stable_set_reduction = true;
                    lp.r[0] = r[0];
                    lp.r[1] = r[1];
                    lp.alpha = alpha;
                }
            }
        } else {
            std::vector<std::vector<bool>> selected(2);
            for (ui i = 0; i < 2; i++) {
                selected[i].resize(n, false);
                for (auto u: vertices[i])
                    selected[i][u] = true;
            }
            std::vector<double> r[2];
            r[0] = lp.r[0];
            r[1] = lp.r[1];
            std::vector<Alpha> alpha = lp.alpha;
            for (ui i = 0; i < lp.edges_count_; i++) {
//            if (!lp.alpha[i].is_selected)
//                continue;
                if (selected[0][alpha[i].id_first] && !selected[1][alpha[i].id_second]) {
                    r[0][alpha[i].id_first] -= 2.0 * sqrt(ratio) * alpha[i].weight_first;
                    r[1][alpha[i].id_second] += 2.0 / sqrt(ratio) * alpha[i].weight_first;
                    alpha[i].weight_first = 0;
                    alpha[i].weight_second = 1;
//                out_degrees[alpha[i].id_first]--;
                } else if (!selected[0][alpha[i].id_first] && selected[1][alpha[i].id_second]) {
                    r[0][alpha[i].id_first] += 2.0 * sqrt(ratio) * alpha[i].weight_second;
                    r[1][alpha[i].id_second] -= 2.0 / sqrt(ratio) * alpha[i].weight_second;
                    alpha[i].weight_first = 1;
                    alpha[i].weight_second = 0;
                    in_degrees[alpha[i].id_second]--;
                } else if (selected[0][alpha[i].id_first] && selected[1][alpha[i].id_second]) {
                    edges.emplace_back(alpha[i].id_first, alpha[i].id_second);
                }
            }
            std::vector<ui> pos(2, 0);
            double min_r = r[cur][tmp_r[cur][0].second];
//        printf("(%f, %f)\n", min_r, rho_c);
            while (cur != best_pos.first || pos[cur] != best_pos.second) {
                if (pos[1] == cnt[1] || (pos[0] != cnt[0] && tmp_r[0][pos[0]] < tmp_r[1][pos[1]])) {
                    cur = 0;
                } else {
                    cur = 1;
                }
//        min_r = std::max(min_r, tmp_r[cur][pos_[cur]].first);
                min_r = std::min(min_r, r[cur][tmp_r[cur][pos[cur]].second]);
                pos[cur]++;
            }
            while (pos[0] < cnt[0] && pos[1] < cnt[1]) {
                if (pos[1] == cnt[1] || (pos[0] != cnt[0] && tmp_r[0][pos[0]] < tmp_r[1][pos[1]])) {
                    cur = 0;
                } else {
                    cur = 1;
                }
//        if (min_r > tmp_r[cur][pos_[cur]].first)
                if (min_r < r[cur][tmp_r[cur][pos[cur]].second])
                    return true;
                pos[cur]++;
            }
//    printf("Stable set.\n");
            VertexID tmp = 1;
            std::vector<std::vector<VertexID>> map(2);
            for (ui i = 0; i < 2; i++) {
                map[i].resize(n, 0);
                for (auto &u: vertices[i])
                    map[i][u] = tmp++;
            }


            ui s = 0, t = tmp;
            FlowNetwork flow = FlowNetwork(tmp + 1);
            for (auto &u: vertices[0]) {
//        flow.addEdge(u, t, rho_c / 2 / sqrt(ratio));
                flow.addEdge(map[0][u], t, rho_c / 2 / sqrt(ratio));
            }
            for (auto &u: vertices[1]) {
//        flow.addEdge(u + n + 1, t, rho_c / 2 * sqrt(ratio));
//        flow.addEdge(s, u + n + 1, in_degrees[u]);
                flow.addEdge(map[1][u], t, rho_c / 2 * sqrt(ratio));
                flow.addEdge(s, map[1][u], in_degrees[u]);
            }
            for (auto &edge: edges) {
//        flow.addEdge(edge.second + n + 1, edge.first + 1, 2.0);
                flow.addEdge(map[1][edge.second], map[0][edge.first], 2.0);
            }
            auto max_flow = flow.getMaxFlow(s, t);
            flag = std::abs(max_flow - edges.size()) > 1e-3;
//        std::vector<bool> is_selected[2];
//        is_selected[0].resize(n, false);
//        is_selected[1].resize(n, false);
//        for (auto &edge: edges) {
//            is_selected[0][edge.first] = true;
//            is_selected[1][edge.second] = true;
//        }
//        vertices[0].clear();
//        vertices[1].clear();
//        for (ui i = 0; i < n; i++) {
//            if (is_selected[0][i])
//                vertices[0].emplace_back(i);
//            if (is_selected[1][i])
//                vertices[1].emplace_back(i);
//        }
//        printf("density: %f, Sc/Tc: %d/%d\n", edges.size() / sqrt((double) vertices[0].size() * vertices[1].size()), vertices[0].size(), vertices[1].size());
//        printf("%e\n", std::abs(max_flow - edges.size()));
//        if (!flag)
//            printf("edges: %d\n", edges.size());
//        lp.r[0] = r[0];
//        lp.r[1] = r[1];
//        lp.alpha = alpha;
            if (edges.size() < m) {
                stable_set_reduction = true;
                lp.r[0] = r[0];
                lp.r[1] = r[1];
                lp.alpha = alpha;
            }
        }
    }
    if (rho > graph.subgraph_density) {
        subgraph.subgraph_density = graph.subgraph_density;
        graph.subgraph_density = rho;
        if (subgraph.is_mapped){
            for (ui i = 0; i < 2; i++)
                for (auto &u: vertices[i])
                    u = subgraph.map[u];
        }
        graph.vertices[0] = vertices[0];
        graph.vertices[1] = vertices[1];
    }
//    if (!flag) {
////        for (ui i = 0; i < 2; i++) {
////            for (auto &u: vertices[i])
////                if (i? !in_degrees[u] : !out_degrees[u])
////                    printf("haha\n");
////        }
//
//        std::vector<bool> is_selected[2];
//        is_selected[0].resize(n, false);
//        is_selected[1].resize(n, false);
//        for (auto &edge: edges) {
//            is_selected[0][edge.first] = true;
//            is_selected[1][edge.second] = true;
//        }
//        vertices[0].clear();
//        vertices[1].clear();
//        for (ui i = 0; i < n; i++) {
//            if (is_selected[0][i])
//                vertices[0].emplace_back(i);
//            if (is_selected[1][i])
//                vertices[1].emplace_back(i);
//        }
//        printf("Sc/Tc: %d/%d\n", vertices[0].size(), vertices[1].size());
//        if (edges.size() / sqrt((double) vertices[0].size() * vertices[1].size()) > graph.subgraph_density){
//            graph.subgraph_density = edges.size() / sqrt((double) vertices[0].size() * vertices[1].size());
//            graph.vertices[0] = vertices[0];
//            graph.vertices[1] = vertices[1];
//        }
//    }

    return flag;
}


bool Verification::directedVWApproVerification(Graph &graph, LinearProgramming &lp,
                                               std::vector<std::vector<VertexID>> &vertices,
                                               double rho, double vw_rho,
                                               double epsilon) {
    if (rho > graph.subgraph_density) {
        graph.subgraph_density = rho;
        graph.vertices[0] = vertices[0];
        graph.vertices[1] = vertices[1];
    }

    auto out_degrees = graph.getOutDegrees();
    auto in_degrees = graph.getInDegrees();
    std::vector<std::pair<VertexID, VertexID>> edges;
    ui n = graph.getVerticesCount();
    double max_r = 0;
    for (VertexID u = 0; u < n; u++) {
        if (out_degrees[u]) {
            max_r = std::max(max_r, lp.r[0][u]);
        }
        if (in_degrees[u]) {
            max_r = std::max(max_r, lp.r[1][u]);
        }
    }

    return vw_rho / max_r < 1 - epsilon / 2;
}


bool Verification::UndirectedlpVerification(Graph &graph, LinearProgramming &lp, FlowNetwork &flow,
                                            std::vector<VertexID> *vertices, bool reduce) {
    std::vector<ui> y;
    std::vector<std::pair<double, VertexID>> tmp;
    ui n = graph.getVerticesCount();
    ui m = graph.getEdgesCount();
    y.resize(n);
    tmp.resize(n);
    for (ui u = 0; u < n; u++) {
        y[u] = 0;
        tmp[u] = std::make_pair(-lp.r[0][u], u);
    }

    for (ui i = 0; i < m; i++) {
        if (lp.r[0][lp.alpha[i].id_first] < lp.r[0][lp.alpha[i].id_second] ||
            (lp.r[0][lp.alpha[i].id_first] == lp.r[0][lp.alpha[i].id_second] &&
             lp.alpha[i].id_first > lp.alpha[i].id_second))
            y[lp.alpha[i].id_first] += 1.0 / lp.weight[lp.alpha[i].id_first];
        else
            y[lp.alpha[i].id_second] += 1.0 / lp.weight[lp.alpha[i].id_second];
    }
    sort(tmp.begin(), tmp.end());
    std::vector<double> val;
    std::vector<ui> sz;
    std::vector<ui> belong;
    belong.resize(n);
    sz.resize(n);
    val.resize(n);
    ui top = -1;
    double eps = 1e-6;
    for (ui u = 0; u < n; u++) {
        top++;
        val[top] = y[tmp[u].second];
        sz[top] = 1;
        while (top >= 1 && (val[top] + val[top - 1]) / (sz[top] + sz[top - 1]) > val[top - 1] / sz[top - 1]) {
            val[top - 1] += val[top];
            sz[top - 1] += sz[top];
            top--;
        }
    }
    ui pos = 0;
    for (ui u = 0; u < n; u++) {
        if (!sz[pos]) pos++;
        sz[pos]--;
        belong[tmp[u].second] = pos;
    }
    for (ui i = 0; i < m; i++) {
        if (belong[lp.alpha[i].id_first] == belong[lp.alpha[i].id_second]) continue;
        if (belong[lp.alpha[i].id_first] < belong[lp.alpha[i].id_second]) {
            lp.r[0][lp.alpha[i].id_second] += lp.alpha[i].weight_first / lp.weight[lp.alpha[i].id_first];
            lp.r[0][lp.alpha[i].id_first] -= lp.alpha[i].weight_first / lp.weight[lp.alpha[i].id_second];
            lp.alpha[i].weight_first = 0;
            lp.alpha[i].weight_second = 1;
        }
        if (belong[lp.alpha[i].id_first] > belong[lp.alpha[i].id_second]) {
            lp.r[0][lp.alpha[i].id_first] += lp.alpha[i].weight_second / lp.weight[lp.alpha[i].id_first];
            lp.r[0][lp.alpha[i].id_second] -= lp.alpha[i].weight_second / lp.weight[lp.alpha[i].id_second];
            lp.alpha[i].weight_second = 0;
            lp.alpha[i].weight_first = 1;
        }
    }
    double minn = 1e20;
    std::vector<ui> rev;
    rev.resize(n);
    pos = 0;
    for (ui u = 0; u < n; u++) {
        rev[tmp[u].second] = u + 1;
        if (belong[tmp[u].second] > pos) {
            pos++;
        }
        if (pos != 0 && lp.r[0][tmp[u].second] >= minn) {
            return true;
        } else if (pos == 0)
            minn = std::min(minn, lp.r[0][tmp[u].second]);
    }

    ui opt_edge_count = 0;
    ui opt_node_count = 0;
    for (ui i = 0; i < vertices[0].size(); i++) opt_node_count += lp.weight[vertices[0][i]];
    for (ui i = 0; i < m; i++) {
        if (belong[lp.alpha[i].id_first] == 0 && belong[lp.alpha[i].id_second] == 0) {
            opt_edge_count++;
        }
    }
    ui S = 0, T = opt_edge_count + vertices[0].size() + 1;
    flow = FlowNetwork(opt_edge_count + vertices[0].size() + 2);
    pos = 0;

    for (ui i = 0; i < m; i++) {
        if (belong[lp.alpha[i].id_first] == 0 && belong[lp.alpha[i].id_second] == 0) {
            pos++;
            flow.addEdge(S, pos, opt_node_count);
            flow.addEdge(pos, opt_edge_count + rev[lp.alpha[i].id_first], opt_node_count);
            flow.addEdge(pos, opt_edge_count + rev[lp.alpha[i].id_second], opt_node_count);
        }
    }
    for (ui i = 1; i <= opt_node_count; i++) {
        flow.addEdge(opt_edge_count + i, T, opt_edge_count * lp.weight[vertices[0][i - 1]]);
    }
    std::cout << "Calling max flow" << std::endl;
    //std::cout<< abs(flow.getMaxFlow(S, T) - 1.0 * opt_edge_count * opt_node_count)<<std::endl;
    double res = abs(flow.getMaxFlow(S, T) - 1.0 * opt_edge_count * opt_node_count);
    if (reduce) {
        std::vector<ui> rev;
        std::vector<ui> id;
        id.resize(n);
        for (ui u = 0; u < n; u++) id[u] = 0;
        ui new_vertices_count_ = 0;
        for (ui i = 0; i < n; i++) {
            if (belong[tmp[i].second] == 0) {
                id[tmp[i].second] = 1 + new_vertices_count_;
                new_vertices_count_++;
            }
        }
        Graph copy(0, new_vertices_count_);
        copy.weight_.resize(new_vertices_count_);
        for (ui i = 0; i < new_vertices_count_; i++) copy.weight_[i] = 1;
        for (ui i = 0; i < n; i++) {
            if (id[i] < 1) continue;
            for (auto v: graph.getNeighbors(i)) {
                if (id[v] < 1 || id[i] < id[v]) continue;
                copy.addUndirectedEdge(id[i] - 1, id[v] - 1);
            }
        }
        std::vector<Alpha> new_alpha;
        for (ui i = 0; i < graph.getEdgesCount(); i++) {
            if (id[lp.alpha[i].id_first] >= 1 && id[lp.alpha[i].id_second] >= 1) {
                lp.alpha[i].id_first = id[lp.alpha[i].id_first] - 1;
                lp.alpha[i].id_second = id[lp.alpha[i].id_second] - 1;
                new_alpha.push_back(lp.alpha[i]);
            }
        }
        lp.alpha = new_alpha;
        if (lp.type_) {
            std::vector<Alpha> new_beta;
            for (ui i = 0; i < graph.getEdgesCount(); i++) {
                if (id[lp.beta[i].id_first] >= 1 && id[lp.beta[i].id_second] >= 1) {
                    lp.beta[i].id_first = id[lp.beta[i].id_first] - 1;
                    lp.beta[i].id_second = id[lp.beta[i].id_second] - 1;
                    new_beta.push_back(lp.beta[i]);
                }
            }
            lp.beta = new_beta;
        }
        lp.nodes_count_ = copy.getVerticesCount();
        lp.edges_count_ = copy.getEdgesCount();
        for (ui i = 0; i < lp.nodes_count_; i++) lp.r[0][i] = 0;
        for (ui i = 0; i < lp.edges_count_; i++) {
            lp.r[0][lp.alpha[i].id_first] += lp.alpha[i].weight_first;
            lp.r[0][lp.alpha[i].id_second] += lp.alpha[i].weight_second;
        }
        graph = copy;
    }
    std::cout << " res = " << res << std::endl;
    std::cout<<1.0 * opt_edge_count / opt_node_count << std::endl;
    return res > eps;
}

bool Verification::UndirectedCoreAppVerification(Graph &graph, CoreApp &ca) {
    std::cout << ca.opt << std::endl;
    ui maxx = 0;
    for(ui i = 0; i < ca.pos; i++){
        if(ca.selected[i]) continue;
        maxx = std::max(maxx, ca.deg[i]);
    }
    if(maxx < ca.k) return false;
    ca.size = ca.pos * 2;
    return true;
}

bool Verification::UndirectedLpAppVerification(Graph &graph, LinearProgramming &lp,
                                               std::vector<VertexID> *vertices, double eps) {
    std::vector<ui> y;
    std::vector<std::pair<double, VertexID>> tmp;
    ui n = graph.getVerticesCount();
    ui m = graph.getEdgesCount();
    y.resize(n);
    tmp.resize(n);
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
    ui pos = 0;
    double sum = 0, opt = 0, weight_sum = 0;
    for (ui i = 0; i < n; i++) {
        sum += y[tmp[i].second];
        weight_sum += lp.weight[tmp[i].second];
        if (sum / weight_sum > opt) {
            opt = sum / weight_sum;
            pos = i + 1;
        }
    }
    double ratio_bound = 0, ratio_real = 0;
    sum = 0, weight_sum = 0;
    for (ui i = 0; i < pos; i++) {
        weight_sum += lp.weight[tmp[i].second];
        sum -= tmp[i].first;
        ratio_bound = std::max(ratio_bound, std::min(sum / (i + 1), 0.5 * i * (i + 1) / weight_sum));
    }
    graph.subgraph_density_lower_bound = opt;
    return 1 + eps < ratio_bound / opt;
}

bool Verification::UndirectedPKMCVerification(Graph &graph, PKMC &pkmc) {
    std::cout << pkmc.hmax << " " << pkmc.s << std::endl;
    bool flag = false;
    ui vertices_count = graph.getVerticesCount();
    PKMC pkmc_new = PKMC();
    pkmc_new.h.resize(vertices_count);
    std::vector<ui> deg = graph.getDegrees();
    std::vector<ui> bucket;
    for (ui i = 0; i < vertices_count; i++) {
        bucket.resize(deg[i] + 1, 0);
        for (ui j = 0; j <= deg[i]; j++) bucket[j] = 0;
        for (auto v: graph.getNeighbors(i)) {
            bucket[std::min(pkmc.h[v], deg[i])] += 1;
        }
        ui res = deg[i];
        while (bucket[res] < res) {
            res--;
            bucket[res] += bucket[res + 1];
        }
        pkmc_new.h[i] = res;
        if (pkmc_new.h[i] < pkmc.h[i]) flag = true;
        pkmc_new.hmax = std::max(pkmc_new.hmax, pkmc_new.h[i]);
    }
    for (ui i = 0; i < vertices_count; i++) {
        pkmc_new.s += (pkmc_new.hmax == pkmc_new.h[i]);
    }
    std::swap(pkmc, pkmc_new);
    if (pkmc.s <= pkmc.hmax) {
        if (flag == false) {
            ui vertex_num = 0, edge_num = 0;
            for (ui i = 0; i < vertices_count; i++) {
                if (pkmc.h[i] == pkmc.hmax) {
                    vertex_num++;
                    for (auto v: graph.getNeighbors(i)) {
                        if (v < i && pkmc.h[v] == pkmc.hmax) edge_num++;
                    }
                }
            }
            std::cout << 1.0 * edge_num / vertex_num << std::endl;
        }
        return flag;
    }
    if (pkmc.hmax == pkmc_new.hmax && pkmc.s == pkmc_new.s) flag = false;
    if (flag == false) {
        ui vertex_num = 0, edge_num = 0;
        for (ui i = 0; i < vertices_count; i++) {
            if (pkmc.h[i] == pkmc.hmax) {
                vertex_num++;
                for (auto v: graph.getNeighbors(i)) {
                    if (v < i && pkmc.h[v] == pkmc.hmax) edge_num++;
                }
            }
        }
        std::cout << 1.0 * edge_num / vertex_num << std::endl;
    }
    return flag;

}

bool Verification::UndirectedFlowAppVerification(Graph &graph, double epsilon){
    epsilon = epsilon / (1 + epsilon);
    double mid = (graph.subgraph_density_lower_bound + graph.subgraph_density_upper_bound) / 2;
    double eps = (mid - graph.subgraph_density_lower_bound) / (2 * mid);
    return eps > epsilon / (3 - 2 * epsilon);
}