//
// Created by yy on 12/1/23.
//

#include "verification.h"

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
    if (edges_count){
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
    if (edges_count){
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
                                     double rho,
                                     double rho_c, double &ratio_o, double &ratio_p, bool &stable_set_reduction,
                                     std::vector<std::pair<VertexID, VertexID>> &edges, double epsilon) {
    double ratio;
//    if (ratios.first < 1 && ratios.second > 1) {
//        ratio = 1;
//    } else if (ratios.second <= 1) {
//        ratio = (ratios.first + ratios.second) / 2;
//    } else if (ratios.first >= 1) {
//        ratio = 2 / (1 / ratios.first + 1 / ratios.second);
//    }
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
//    if (!cnt[0] || !cnt[1])
//        return true;
    for (ui i = 0; i < 2; i++)
        sort(tmp_r[i].begin(), tmp_r[i].end());

    ui cur = tmp_r[0][0] > tmp_r[1][0] ? 1 : 0;
    bool flag;
//    if (-tmp_r[cur][0].first < graph.subgraph_density * sqrt(1 + epsilon)) {
//        double t = -tmp_r[cur][0].first / graph.subgraph_density / sqrt(1 + epsilon);
//        ratio_o = (2 * ratio - ratio * t * t - 2 * sqrt(ratio * ratio - ratio * ratio * t * t)) / (t * t);
//        ratio_p = ratio * ratio / ratio_o;
//        if (ratio_o > ratio_p)
//            std::swap(ratio_o, ratio_p);
//        flag = false;
//    }
    if (epsilon) {
//        flag = -tmp_r[cur][0].first / rho_c > sqrt(1 + epsilon);
//        flag = -tmp_r[cur][0].first > sqrt(1 + epsilon) * rho_c;
//        ratio_o = std::min(ratio_o, ratio / (1 + epsilon));
//        ratio_p = std::max(ratio_p, ratio * (1 + epsilon));
        if (-tmp_r[cur][0].first / rho_c <= sqrt(1 + epsilon)) {
            flag = false;
            ratio_o = std::min(ratio_o, ratio / (1 + epsilon));
            ratio_p = std::max(ratio_p, ratio * (1 + epsilon));
        } else if (-tmp_r[cur][0].first / rho_c <= 1 + epsilon && ratio_o < ratio / (1 + epsilon) && ratio_p > ratio * (1 + epsilon)) {
            flag = false;
        } else
            flag = true;
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
        double min_r = -tmp_r[cur][0].first;
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
        if (edges.size() < m) {
            stable_set_reduction = true;
            lp.r[0] = r[0];
            lp.r[1] = r[1];
            lp.alpha = alpha;
        }
    }
    if (rho > graph.subgraph_density) {
        graph.subgraph_density = rho;
        graph.vertices[0] = vertices[0];
        graph.vertices[1] = vertices[1];
    }

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
                                            std::vector<VertexID> *vertices) {
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

    for(ui i = 0; i < m; i++){
        if(lp.r[0][lp.alpha[i].id_first] < lp.r[0][lp.alpha[i].id_second] || 
            (lp.r[0][lp.alpha[i].id_first] == lp.r[0][lp.alpha[i].id_second] && lp.alpha[i].id_first > lp.alpha[i].id_second))
            y[lp.alpha[i].id_first]++;
        else
            y[lp.alpha[i].id_second]++;
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
    for(ui i = 0; i < m; i++){
        if(belong[lp.alpha[i].id_first] == belong[lp.alpha[i].id_second]) continue;
        if(belong[lp.alpha[i].id_first] < belong[lp.alpha[i].id_second]){
            lp.r[0][lp.alpha[i].id_second] += lp.alpha[i].weight_first;
            lp.r[0][lp.alpha[i].id_first] -= lp.alpha[i].weight_first;
            lp.alpha[i].weight_first = 0;
            lp.alpha[i].weight_second = 1;
        }
        if(belong[lp.alpha[i].id_first] > belong[lp.alpha[i].id_second]){
            lp.r[0][lp.alpha[i].id_first] += lp.alpha[i].weight_second;
            lp.r[0][lp.alpha[i].id_second] -= lp.alpha[i].weight_second;
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
    ui opt_node_count = vertices[0].size();

    for(ui i = 0; i < m; i++){
        if(belong[lp.alpha[i].id_first] == 0 && belong[lp.alpha[i].id_second] == 0){
            opt_edge_count++;
        }
    }
    ui S = 0, T = opt_edge_count + opt_node_count + 1;
    flow = FlowNetwork(opt_edge_count + opt_node_count + 2);
    pos = 0;

    for(ui i = 0; i < m; i++){
        if(belong[lp.alpha[i].id_first] == 0 && belong[lp.alpha[i].id_second] == 0){
            pos++;
            flow.addEdge(S,pos,opt_node_count);
            flow.addEdge(pos, opt_edge_count + rev[lp.alpha[i].id_first], opt_node_count);
            flow.addEdge(pos, opt_edge_count + rev[lp.alpha[i].id_second], opt_node_count);
        }
    }
    for (ui i = 1; i <= opt_node_count; i++) {
        flow.addEdge(opt_edge_count + i, T, opt_edge_count);
    }
    return abs(flow.getMaxFlow(S, T) - 1.0 * opt_edge_count * opt_node_count) > eps;
}
//bool Verification::UndirectedCoreAppVerification(Graph &graph, CoreApp &ca){
//    std::vector<ui> deg = graph.getDegrees();
//    if(ca.pos + 1 == ca.nodes_count || deg[ca.id[ca.pos + 1]] < ca.k) return false;
//    ca.pos = 2 * ca.pos + 1;
//    if(ca.pos + 1 > ca.nodes_count){
//        ca.pos = ca.nodes_count - 1;
//    }
//    while(ca.pos + 1 < ca.nodes_count && deg[ca.id[ca.pos]] == deg[ca.id[ca.pos + 1]]) ca.pos++;
//    return true;
//}