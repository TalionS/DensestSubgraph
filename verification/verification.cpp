//
// Created by yy on 12/1/23.
//

#include "verification.h"

bool Verification::flowExactVerification(Graph &graph, double l, double r) {
    ui n = graph.getVerticesCount();
    double bias = 1.0 / sqrt((double) n * (n - 1)) - 1.0 / n;

    if (r - l > bias) return true;
    else return false;
}

bool Verification::UndirectedflowExactVerification(Graph &graph, double l, double r) {
    ui n = graph.getVerticesCount();
    double bias = 1.0 / ((double) n * (n - 1));
    if (r - l > bias) return true;
    else return false;
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
        if(lp.r[0][lp.alpha[0][i].id_first] < lp.r[0][lp.alpha[0][i].id_second] || 
            (lp.r[0][lp.alpha[0][i].id_first] == lp.r[0][lp.alpha[0][i].id_second] && lp.alpha[0][i].id_first > lp.alpha[0][i].id_second))
            y[lp.alpha[0][i].id_first]++;
        else
            y[lp.alpha[0][i].id_second]++;
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
        if(belong[lp.alpha[0][i].id_first] == belong[lp.alpha[0][i].id_second]) continue;
        if(belong[lp.alpha[0][i].id_first] < belong[lp.alpha[0][i].id_second]){
            lp.r[0][lp.alpha[0][i].id_second] += lp.alpha[0][i].weight_first;
            lp.r[0][lp.alpha[0][i].id_first] -= lp.alpha[0][i].weight_first;
            lp.alpha[0][i].weight_first = 0;
            lp.alpha[0][i].weight_second = 1;
        }
        if(belong[lp.alpha[0][i].id_first] > belong[lp.alpha[0][i].id_second]){
            lp.r[0][lp.alpha[0][i].id_first] += lp.alpha[0][i].weight_second;
            lp.r[0][lp.alpha[0][i].id_second] -= lp.alpha[0][i].weight_second;
            lp.alpha[0][i].weight_second = 0;
            lp.alpha[0][i].weight_first = 1;
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
        if(belong[lp.alpha[0][i].id_first] == 0 && belong[lp.alpha[0][i].id_second] == 0){
            opt_edge_count++;
        }
    }
    ui S = 0, T = opt_edge_count + opt_node_count + 1;
    flow = FlowNetwork(opt_edge_count + opt_node_count + 2);
    pos = 0;
    for(ui i = 0; i < m; i++){
        if(belong[lp.alpha[0][i].id_first] == 0 && belong[lp.alpha[0][i].id_second] == 0){
            pos++;
            flow.addEdge(S,pos,opt_node_count);
            flow.addEdge(pos, opt_edge_count + rev[lp.alpha[0][i].id_first], opt_node_count);
            flow.addEdge(pos, opt_edge_count + rev[lp.alpha[0][i].id_second], opt_node_count);
        }
    }
    for(ui i = 1; i <= opt_node_count; i++){
        flow.addEdge(opt_edge_count + i, T, opt_edge_count);
    }
    return abs(flow.getMaxFlow(S,T) - 1.0 * opt_edge_count * opt_node_count) > eps;
}