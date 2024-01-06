//
// Created by yy on 12/1/23.
//

#include "verification.h"

bool Verification::flowExactVerification(Graph &graph, XYCore xy_core, double l, double r) {
    ui n = graph.getVerticesCount();
    double bias = 1.0 / sqrt((double) n * (n - 1)) - 1.0 / n;

//    double bias = 1 / sqrt((double) (n - xy_core.bin[0][1]) * (n - xy_core.bin[1][1]) - std::min(n - xy_core.bin[0][1], n - xy_core.bin[1][1]))
//    - 1 / sqrt((double) (n - xy_core.bin[0][1]) * (n - xy_core.bin[1][1]));
    if (r - l > bias) return true;
    else return false;
}

bool Verification::UndirectedflowExactVerification(Graph &graph, double l, double r) {
    ui n = graph.getVerticesCount();
    double bias = 1.0 / ((double) n * (n - 1));
    if (r - l > bias) return true;
    else return false;
}

bool Verification::directedLPExactVerification(Graph &graph,
                                               Graph &x_y_core,
                                               LinearProgramming &lp,
                                               std::pair <ui, ui> best_pos,
                                               std::vector<std::vector<VertexID>> &vertices,
                                               std::pair<double, double> ratios,
                                               double rho_c) {
    double ratio = (ratios.first + ratios.second) / 2;
//    double ratio = ratios.first / ratios.second;
    std::vector<std::vector<std::pair<double, VertexID>>> tmp_r(2);
    std::vector<ui> cnt(2, 0);
    auto out_degrees = x_y_core.getOutDegrees();
    auto in_degrees = x_y_core.getInDegrees();
    std::vector<std::pair<VertexID, VertexID>> edges;
    ui n = x_y_core.getVerticesCount();
    ui m = x_y_core.getEdgesCount();
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
    for (ui i = 0; i < 2; i++)
        sort(tmp_r[i].begin(), tmp_r[i].end());
    std::vector<std::vector<bool>> selected(2);
    for(ui i = 0; i < 2; i++) {
        selected[i].resize(n, false);
        for (auto u: vertices[i])
            selected[i][u] = true;
    }
    for (ui i = 0; i < m; i++){
        if (selected[0][lp.alpha[i].id_first] && !selected[1][lp.alpha[i].id_second]){
            lp.r[0][lp.alpha[i].id_first] -= 2.0 * sqrt(ratio) * lp.alpha[i].weight_first;
            lp.r[1][lp.alpha[i].id_second] += 2.0 / sqrt(ratio) * lp.alpha[i].weight_first;
            lp.alpha[i].weight_first = 0;
            lp.alpha[i].weight_second = 1;
        }
        else if(!selected[0][lp.alpha[i].id_first] && selected[1][lp.alpha[i].id_second]){
            lp.r[0][lp.alpha[i].id_first] += 2.0 * sqrt(ratio) * lp.alpha[i].weight_second;
            lp.r[1][lp.alpha[i].id_second] -= 2.0 / sqrt(ratio) * lp.alpha[i].weight_second;
            lp.alpha[i].weight_first = 1;
            lp.alpha[i].weight_second = 0;
            in_degrees[lp.alpha[i].id_second]--;
        }
        else if(selected[0][lp.alpha[i].id_first] && selected[1][lp.alpha[i].id_second])
            edges.emplace_back(std::make_pair(lp.alpha[i].id_first, lp.alpha[i].id_second));
    }
    std::vector<ui> pos(2, 0);
    if (!cnt[0] || !cnt[1])
        return false;
    ui cur = tmp_r[0][0] > tmp_r[1][0] ? 1 : 0;
    double min_r = tmp_r[cur][0].first;
    while(cur != best_pos.first || pos[cur] != best_pos.second){
        if (pos[1] == cnt[1] || tmp_r[0][pos[0]] < tmp_r[1][pos[1]]){
            cur = 0;
        }
        else{
            cur = 1;
        }
        min_r = std::max(min_r, tmp_r[cur][pos[cur]].first);
        pos[cur]++;
    }
    while(pos[0] < cnt[0] && pos[1] < cnt[1]){
        if (pos[1] == cnt[1] || tmp_r[0][pos[0]] < tmp_r[1][pos[1]]){
            cur = 0;
        }
        else{
            cur = 1;
        }
        if (min_r > tmp_r[cur][pos[cur]].first)
            return true;
        pos[cur]++;
    }
//    printf("Stable set.\n");
    ui s = 0, t = 2 * n - 1;
    FlowNetwork flow = FlowNetwork(2 * n + 2);
    for (auto u: vertices[0]){
        flow.addEdge(u, t, rho_c / 2 / sqrt(ratio));
    }
    for (auto u: vertices[1]){
        flow.addEdge(u + n + 1, t, sqrt(ratio) * rho_c / 2);
        flow.addEdge(s, u + n + 1, in_degrees[u]);
    }
    for (auto edge: edges){
        flow.addEdge(edge.second + n + 1, edge.first + 1, 2);
    }
    auto max_flow = flow.getMaxFlow(s, t);
    bool flag = std::abs(max_flow - edges.size()) > 1;
    if (edges.size() / sqrt(vertices[0].size() * vertices[1].size()) > graph.subgraph_density){
        graph.subgraph_density = edges.size() / sqrt(vertices[0].size() * vertices[1].size());
        graph.vertices = &vertices[0];
    }
    if(edges.size() < m){
        Graph stable_subgraph = Graph(true, n);
        for (auto edge : edges)
            stable_subgraph.addDirectedEdge(edge.first, edge.second);
        x_y_core = stable_subgraph;
    }
    return flag;
}

bool Verification::UndirectedlpVerification(Graph &graph, LinearProgramming &lp, FlowNetwork &flow, std::vector<VertexID> *vertices) {
    std::vector<ui> y;
    std::vector<std::pair<double,VertexID>> tmp;
    ui n = graph.getVerticesCount();
    ui m = graph.getEdgesCount();
    y.resize(n);
    tmp.resize(n);
    for(ui u = 0; u < n; u++){
        y[u] = 0;
        tmp[u] = std::make_pair(-lp.r[0][u],u);
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
    for(ui u = 0; u < n; u++){
        top++;
        val[top] = y[tmp[u].second];
        sz[top] = 1;
        while(top >= 1 && (val[top] + val[top - 1]) / (sz[top] + sz[top - 1]) > val[top - 1] / sz[top - 1]){
            val[top - 1] += val[top];
            sz[top - 1] += sz[top];
            top--;
        }
    }
    ui pos = 0;
    for(ui u = 0; u < n; u++){
        if(!sz[pos]) pos++;
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
    for(ui u = 0; u < n; u++){
        rev[tmp[u].second] = u + 1;
        if(belong[tmp[u].second] > pos){
            pos++;
        }
        if(pos != 0 && lp.r[0][tmp[u].second] >= minn){
            return true;
        }
        else if(pos == 0)
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
    for(ui i = 1; i <= opt_node_count; i++){
        flow.addEdge(opt_edge_count + i, T, opt_edge_count);
    }
    return abs(flow.getMaxFlow(S,T) - 1.0 * opt_edge_count * opt_node_count) > eps;
}