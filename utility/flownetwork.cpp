//
// Created by yy on 11/19/23.
//

#include "flownetwork.h"
#include "graph.h"
#include <iostream>
#include <queue>

FlowEdge::FlowEdge(VertexID from, VertexID to, double capacity, double flow, VertexID index) : from(from), to(to), index(index), capacity(capacity), flow(flow) {};

//FlowNetwork::FlowNetwork() : nodes_count_(0){
//
//}

FlowNetwork::FlowNetwork(ui vertices_count) : nodes_count_(vertices_count){
    adj_.clear();
    adj_.resize(vertices_count);
//    vertices.resize(vertices_count);
//    for (int i = 0; i < vertices_count; ++i) {
//        vertices[i] = i;
//    }
}


void FlowNetwork::addEdge(VertexID from, VertexID to, double capacity) {
//    printf("%d %d %f\n", from, to, capacity);
    FlowEdge forward_edge(from, to, capacity, 0, adj_[to].size());
    adj_[from].push_back(forward_edge);
    if(from == to){
        adj_[from].back().index++;
    }
    FlowEdge reverse_edge(to, from, 0, 0, adj_[from].size() - 1);
    adj_[to].push_back(reverse_edge);
}

void FlowNetwork::enqueue(VertexID v) {
    if(!active_[v] && excess_[v] > 0 && dist_[v] < nodes_count_){
        active_[v] = true;
        height_bucket_[dist_[v]].push_back(v);
        highest_active_height_ = std::max(highest_active_height_, dist_[v]);
    }
}

void FlowNetwork::push(FlowEdge &e) {
    double flow_to_push = std::min(excess_[e.from], e.capacity - e.flow);
    if(dist_[e.from] == dist_[e.to] + 1 && flow_to_push > 0){
        e.flow += flow_to_push;
        adj_[e.to][e.index].flow -= flow_to_push;
        excess_[e.to] += flow_to_push;
        excess_[e.from] -= flow_to_push;
        enqueue(e.to);
    }
}

void FlowNetwork::gap(VertexID k) {
    for(VertexID v = 0; v < nodes_count_; v++){
        if(dist_[v] >= k){
            count_[dist_[v]]--;
            dist_[v] = std::max(dist_[v], nodes_count_);
            count_[dist_[v]]++;
            enqueue(v);
        }
    }
}

void FlowNetwork::relabel(VertexID v) {
    count_[dist_[v]]--;
    dist_[v] = nodes_count_;
    for(auto e : adj_[v]){
        if(e.capacity - e.flow > 0){
            dist_[v] = std::min(dist_[v], dist_[e.to] + 1);
        }
    }
    count_[dist_[v]]++;
    enqueue(v);
}

void FlowNetwork::discharge(VertexID v) {
    for (auto &e : adj_[v]) {
        if (excess_[v] > 0) {
            push(e);
        } else {
            break;
        }
    }

    if (excess_[v] > 0) {
        if (count_[dist_[v]] == 1) {
            gap(dist_[v]);
        } else {
            relabel(v);
        }
    }
}

double FlowNetwork::getMaxFlow(VertexID s, VertexID t) {
    dist_ = std::vector<ui>(nodes_count_, 0);
    excess_ = std::vector<double>(nodes_count_, 0);
    count_ = std::vector<ui>(nodes_count_ + 1, 0);
    active_ = std::vector<bool>(nodes_count_, false);
    height_bucket_ = std::vector<std::vector<ui>>(nodes_count_);
    highest_active_height_ = 0;

    for (auto &e: adj_[s]) {
        excess_[s] += e.capacity;
    }

    count_[0] = nodes_count_;
    enqueue(s);
    active_[t] = true;

    while (highest_active_height_ >= 0) {
        if (!height_bucket_[highest_active_height_].empty()) {
            ui v = height_bucket_[highest_active_height_].back();
            height_bucket_[highest_active_height_].pop_back();
            active_[v] = false;
            discharge(v);
        } else {
            if (!highest_active_height_)
                break;
            highest_active_height_--;
        }
    }
//    for(auto &dis: dist_)
//        printf("%d ", dis);
//    printf("\n");
    return excess_[t];
}

void FlowNetwork::getMinCut(VertexID s, VertexID t, std::vector<VertexID> &S) {
//    S.clear();
//    T.clear();
//    std::vector<bool> visited(nodes_count_, false);
//
//    std::queue<VertexID> q;
//    q.push(s);
//    visited[s] = true;
//    while (!q.empty()) {
//        VertexID u = q.front();
//        q.pop();
//
//        for (const auto &e : adj_[u]) {
////            printf("(%d, %d, %f, %f)\n", e.from, e.to, e.flow, e.capacity);
//            if (!visited[e.to] && e.capacity - e.flow > 0) {
//                visited[e.to] = true;
//                q.push(e.to);
//            }
//        }
//    }
//
//    for (VertexID i = 0; i < nodes_count_; ++i) {
//        if (visited[i]) {
//            S.push_back(i);
//        } else {
//            T.push_back(i);
//        }
//    }
    S.clear();
    for (int v = 0; v < nodes_count_; v++) {
//        if (dist_[v] >= nodes_count_) {
        if (dist_[v] >= nodes_count_) {
            S.push_back(v);
        }
    }
}

Edge::Edge(ui from, ui to, double flow, ui rev = 0): from(from), to(to), flow(flow), rev(rev){};
void Dinic::Init(ui n){
    s = 0, t = n - 1;
    e.resize(n);
    dis.resize(n);
    for(ui i = 0; i < n; i++) e[i].clear();
}
void Dinic::addEdge(ui from, ui to, double flow){
    e[from].push_back(Edge(from, to, flow));
    e[to].push_back(Edge(to,from,0));
    e[from].back().rev = e[to].size() - 1;
    e[to].back().rev = e[from].size() - 1;
}
bool Dinic::bfs(){
    for(ui i = s; i <= t; i++) dis[i] = -1;
    std::queue<ui> q;
    q.push(s);
    dis[s] = 0;
    while(q.size()){
        ui u = q.front();q.pop();
        for(ui i = 0; i < e[u].size(); i++){
            if(e[u][i].flow <= 0) continue;
            ui v = e[u][i].to;
            if(dis[v] == -1) dis[v] = dis[u] + 1, q.push(v);
        }
    } 
    return dis[t] != -1;
}

double Dinic::dfs(ui u, double flow){
    if(u == t) return flow;
    double res = flow;
    for(ui i = 0; i < e[u].size(); i++){
        if(e[u][i].flow <= 0) continue;
        ui v = e[u][i].to;
        if(dis[v] == dis[u] + 1){
            double k = dfs(v, std::min(flow, e[u][i].flow));
            if(k <= 0) dis[v] = -1;
            else{
                flow -= k;
                e[u][i].flow -= k;
                e[v][e[u][i].rev].flow += k;
            }
            if(!flow) return res; 
        }
    }
    return res - flow;
}

double Dinic::solve(double iter){
    double h = 0;
    double res = 0;
    while(bfs() && h < iter){
        h++;
        res += dfs(0, 1e10);
    }
    return res;
}