//
// Created by yy on 23-11-13.
//

#include "graph.h"
#include "heap.h"
#include "app.h"
#include <fstream>
#include <iostream>
#include <assert.h>

Graph::Graph(bool is_directed, ui n, bool is_weighted) :
        is_directed_(is_directed),
        is_weighted_(is_weighted),
        vertices_count_(n),
        edges_count_(0),
        subgraph_density(0),
        subgraph_density_lower_bound(0),
        subgraph_density_upper_bound(0) {
    if (is_directed_) {
        adj_ = new std::vector<std::vector<VertexID>>[2];
        deg_ = new std::vector<ui>[2];
        vertices = new std::vector<VertexID>[2];
        for (int i = 0; i < 2; i++) {
            adj_[i].resize(static_cast<unsigned long>(n));
            deg_[i].resize(static_cast<unsigned long>(n));
//            vertices[i].resize(static_cast<unsigned long>(n));
        }
//        vertex_ids.resize(2);
//        map.clear();
//        for (ui i = 0; i < n; i++)
//            map.emplace_back(i);
    } else {
        adj_ = new std::vector<std::vector<VertexID>>[1];
        deg_ = new std::vector<ui>[1];
        vertices = new std::vector<VertexID>[1];
        adj_[0].resize(static_cast<unsigned long>(n));
        weight_.resize(static_cast<double>(n));
        deg_[0].resize(static_cast<unsigned long>(n));
        vertices[0].resize(static_cast<unsigned long>(n));
    }
}

Graph &Graph::operator=(const Graph &other) {
    if (this != &other) {  // Avoid self-assignment
        // Copy data members
        is_directed_ = other.is_directed_;
        is_weighted_ = other.is_weighted_;
        vertices_count_ = other.vertices_count_;
        edges_count_ = other.edges_count_;
        subgraph_density = other.subgraph_density;
        subgraph_density_lower_bound = other.subgraph_density_lower_bound;
        subgraph_density_upper_bound = other.subgraph_density_upper_bound;
        weight_ = other.weight_;
        map = other.map;
        is_mapped = other.is_mapped;
//        vertex_ids = other.vertex_ids;

        // Resize and copy adjacency lists, degree vectors, and vertices vectors
        if (is_directed_) {
            adj_[0] = other.adj_[0];
            adj_[1] = other.adj_[1];
            deg_[0] = other.deg_[0];
            deg_[1] = other.deg_[1];
            vertices[0] = other.vertices[0];
            vertices[1] = other.vertices[1];
//            map = other.map;
        } else {
            adj_[0] = other.adj_[0];
            deg_[0] = other.deg_[0];
            vertices[0] = other.vertices[0];
        }
    }
    return *this;
}

Graph::~Graph() {
    delete[] adj_;
    delete[] deg_;
    delete[] vertices;
}

void Graph::init() {
    if (is_directed_) {
        ui max_out = 0;
        VertexID pos_max_out, pos_max_in;
        ui max_in = 0;
        for (VertexID u = 0; u < vertices_count_; u++) {
            if (deg_[0][u] > max_out) {
                max_out = deg_[0][u];
                pos_max_out = u;
            }
            if (deg_[1][u] > max_in) {
                max_in = deg_[1][u];
                pos_max_in = u;
            }
        }
        subgraph_density_upper_bound = std::max(max_out, max_in);
        subgraph_density = sqrt(subgraph_density_upper_bound);
        ui cur = max_out > max_in? 0 : 1;
        vertices[cur].emplace_back(cur? pos_max_in: pos_max_out);
        vertices[1 - cur] = adj_[cur][cur? pos_max_in: pos_max_out];
    }
}
void Graph::loadGraphFromFile(const std::string &dir, bool is_sample, double sample_rate) {
    std::ifstream file(dir);
    if (!file.is_open()) {
        std::cout << "Unable to open file " << dir << "." << std::endl;
        exit(-1);
    }

    ui edges_count;
    file >> vertices_count_ >> edges_count;
    std::srand(42);
    std::vector<ui> id;
    std::vector<ui> in_sample;
    if(is_sample){
        id.resize(vertices_count_);
        in_sample.resize(vertices_count_);
        ui sample_num = vertices_count_ * sample_rate;
        for(ui i = 0; i < vertices_count_; i++) id[i] = i, in_sample[i] = 0;
        std::random_shuffle(id.begin(), id.end());
        for(ui i = 0; i < sample_num; i++) in_sample[id[i]] = 1;
    }
    

    if (is_directed_) {
        for (int i = 0; i < 2; i++) {
            adj_[i].resize(static_cast<VertexID>(vertices_count_));
            deg_[i].resize(static_cast<VertexID>(vertices_count_), 0);
        }
        VertexID begin, end;
        while (file >> begin >> end) {
            if(is_sample && (!in_sample[begin] || !in_sample[end])) continue;
            addDirectedEdge(begin, end);
        }
//        map.resize(static_cast<VertexID>(vertices_count_));
//        for (ui u = 0; u < vertices_count_; u++)
//            map[u] = u;
    } else {
        adj_[0].resize(static_cast<VertexID>(vertices_count_));
        deg_[0].resize(static_cast<VertexID>(vertices_count_), 0);
        VertexID begin, end;
        while (file >> begin >> end) {
            if(is_sample && (!in_sample[begin] || !in_sample[end])) continue;
            addUndirectedEdge(begin, end);
        }
    }

//    assert(edges_count_ == edges_count);
    file.close();
    if (is_directed_) {
        ui max_out = *std::max_element(deg_[0].begin(), deg_[0].end());
        ui max_in = *std::max_element(deg_[1].begin(), deg_[1].end());
        subgraph_density_upper_bound = std::max(max_out, max_in);
        subgraph_density = sqrt(subgraph_density_upper_bound);
    } else {
        subgraph_density_upper_bound = *std::max_element(deg_[0].begin(), deg_[0].end());
    }
    is_mapped = false;
};

void Graph::sample(double p){

}
void Graph::removeMultiEdges(Graph &graph) {
    Graph copy(is_directed_, vertices_count_);
    if (is_directed_) {
        for (ui u = 0; u < vertices_count_; u++) {
            sort(adj_[0][u].begin(), adj_[0][u].end());
            if (!adj_[0][u].empty()) {
                auto pre_v = graph.getVerticesCount();
                for (auto v: adj_[0][u]) {
                    if (v != pre_v) {
                        copy.addDirectedEdge(u, v);
                    }
                    pre_v = v;
                }
            }
        }
        ui max_out = *std::max_element(copy.deg_[0].begin(), copy.deg_[0].end());
        ui max_in = *std::max_element(copy.deg_[1].begin(), copy.deg_[1].end());
        copy.subgraph_density_upper_bound = std::max(max_out, max_in);
        copy.subgraph_density = sqrt(subgraph_density_upper_bound);
//        copy.vertices[0].clear();
//        copy.vertices[1].clear();
//        for (VertexID u = 0; u < vertices_count_; u++) {
//            if (copy.deg_[0][u])
//                copy.vertex_ids[0].emplace_back(u);
//            if (copy.deg_[1][u])
//                copy.vertex_ids[1].emplace_back(u);
//        }
    } else {
        for (ui u = 0; u < vertices_count_; u++) {
            sort(adj_[0][u].begin(), adj_[0][u].end());
            if (!adj_[0][u].empty()) {
                auto pre_v = graph.getVerticesCount();
                for (auto v: adj_[0][u]) {
                    if (v >= u)
                        break;
                    if (v != pre_v) {
                        copy.addUndirectedEdge(u, v);
                    }
                    pre_v = v;
                }
            }
        }
    }
    graph = copy;
    if(!is_weighted_){
//        for(ui i = 0; i < vertices_count_; i++){
//            weight_[i] = 1;
//        }
//
        weight_.assign(vertices_count_, 1);
        //weight_[3] = 3;
    }
}

void Graph::coreOrder(Graph &graph){
    Graph copy(is_directed_, vertices_count_);
    std::vector<ui> *bucket = new std::vector<ui>[vertices_count_];
    for(ui i = 0 ; i < vertices_count_; i++) bucket[i].clear();
    std::vector<ui> deg = graph.getDegrees();
    std::vector<ui> id;
    id.resize(vertices_count_);
    for(ui i = 0; i < vertices_count_; i++){
        bucket[deg[i]].push_back(i);
    }
    ui pos = 0;
    for(int i = vertices_count_ - 1; i >= 0; i--){
        for(auto u : bucket[i]){
            id[u] = pos++;
        }
    }
    std::vector<ui> *to = new std::vector<ui>[vertices_count_];
    for(ui i = 0 ; i < vertices_count_; i++) to[i].clear();
    for(ui i = 0; i < vertices_count_; i++){
        for(auto v : graph.getNeighbors(i)){
            to[id[v]].push_back(id[i]);
        }
    }
    copy.deg_[0].resize(vertices_count_);
    copy.adj_[0].resize(vertices_count_);
    copy.weight_.resize(vertices_count_);
    for(ui i = 0; i < vertices_count_; i++) copy.weight_[i] = 1;
    for(ui i = 0; i < vertices_count_; i++){
        for(auto v: to[i]){
            copy.edges_count_++;
            copy.deg_[0][v]++;
            copy.adj_[0][v].push_back(i);
        }
    }
    copy.edges_count_ /= 2;
    graph = copy;
}

std::vector<ui> Graph::CoreOrder(Graph &graph){
    std::vector<ui> *bucket = new std::vector<ui>[vertices_count_];
    for(ui i = 0 ; i < vertices_count_; i++) bucket[i].clear();
    std::vector<ui> deg = graph.getDegrees();
    std::vector<ui> id;
    id.resize(vertices_count_);
    for(ui i = 0; i < vertices_count_; i++){
        bucket[deg[i]].push_back(i);
        id[i] = vertices_count_;
    }
    ui pos = 0;
    for(ui i = 0; i < vertices_count_; i++){
        for(ui j = 0; j < bucket[i].size(); j++){
            ui u = bucket[i][j];
            if(id[u] != vertices_count_) continue;
            deg[u] = i;
            id[u] = pos++;
            for(auto v: graph.getNeighbors(u)){
                if(id[v] != vertices_count_) continue;
                bucket[--deg[v]].push_back(v);
            }
        }
    }
    return deg;
}

void Graph::coreReduce(Graph &graph, ui k, bool weighted){
    if(!weighted){
        std::vector<ui> *bucket = new std::vector<ui>[vertices_count_];
        for(ui i = 0 ; i < vertices_count_; i++) bucket[i].clear();
        std::vector<ui> deg = graph.getDegrees();
        std::vector<ui> id;
        id.resize(vertices_count_);
        for(ui i = 0; i < vertices_count_; i++){
            bucket[deg[i]].push_back(i);
            id[i] = vertices_count_;
        }
        double low = 0;
        for(ui i = 0; i < vertices_count_; i++){
            for(ui j = 0; j < bucket[i].size(); j++){
                ui u = bucket[i][j];
                if(id[u] != vertices_count_) continue;
                low = std::max(low, (double) i / 2);
                k = std::max(k, i / 2);
                id[u] = i;
                for(auto v : graph.getNeighbors(u)){
                    if(id[v] != vertices_count_) continue;
                    bucket[--deg[v]].push_back(v);
                }
            }
        }
        ui new_vertices_count_ = 0;
        for(ui i = 0; i < vertices_count_; i++){
            if(id[i] >= k){
                id[i] = k + new_vertices_count_;
                new_vertices_count_++;
            }
        }
        Graph copy(is_directed_, new_vertices_count_);
        copy.deg_[0].resize(new_vertices_count_);
        copy.adj_[0].resize(new_vertices_count_);
        copy.weight_.resize(new_vertices_count_);
        copy.subgraph_density_lower_bound = std::max(low, graph.subgraph_density_lower_bound);
        for(ui i = 0; i < new_vertices_count_; i++) copy.weight_[i] = 1;
        for(ui i = 0; i < vertices_count_; i++){
            if(id[i] < k) continue;
            for(auto v: graph.getNeighbors(i)){
                if(id[v] < k || id[i] < id[v]) continue;
                copy.addUndirectedEdge(id[i] - k, id[v] - k);
            }
        }
        deg = copy.getDegrees();
        double up = 0;
        for(ui i = 0; i < new_vertices_count_; i++){
            up = std::max(up, (double) deg[i]);
        }
        copy.subgraph_density_upper_bound = up;
        graph = copy;
    }
    else{
        BinaryHeap heap;
        heap.resize(vertices_count_);
        std::vector<ui> deg = graph.getDegrees();
        std::vector<double> h;
        h.resize(vertices_count_);
        std::vector<ui> id;
        id.resize(vertices_count_);
        for(ui i = 0; i < vertices_count_; i++) h[i] = deg[i] / weight_[i], id[i] = 0;
        heap.init(h);
        graph.subgraph_density_lower_bound = 0;
        double low = graph.subgraph_density_lower_bound;
        for(ui i = 0; i < vertices_count_; i++){
            HeapElement tmp = heap.pop();
            h[tmp.id] = tmp.val;
            id[tmp.id] = 1;
            low = std::max(low, tmp.val / 2);
            for(auto v: graph.getNeighbors(tmp.id)){
                if(id[v] != 0) continue;
                heap.dec(v, 1.0 / weight_[v], tmp.val);
            }
        }
        
        ui new_vertices_count_ = 0;
        for(ui i = 0; i < vertices_count_; i++){
            if(h[i] >= low) id[i] = new_vertices_count_++;
            else id[i] = vertices_count_ + 1;
        }
        Graph copy(is_directed_, new_vertices_count_);
        copy.deg_[0].resize(new_vertices_count_);
        copy.adj_[0].resize(new_vertices_count_);
        copy.weight_.resize(new_vertices_count_);
        copy.map.resize(new_vertices_count_);
        copy.subgraph_density_lower_bound = std::max(low, graph.subgraph_density_lower_bound);
        new_vertices_count_ = 0;
        for(ui i = 0; i < vertices_count_; i++){
            if(id[i] != vertices_count_ + 1){
                copy.map[new_vertices_count_] = i;
                copy.weight_[new_vertices_count_++] = weight_[i];
            }
        }
        for(ui i = 0; i < vertices_count_; i++){
            if(id[i] == vertices_count_ + 1) continue;
            for(auto v: graph.getNeighbors(i)){
                if(id[v] ==  vertices_count_ + 1|| id[i] < id[v]) continue;
                copy.addUndirectedEdge(id[i], id[v]);
            }
        }
        deg = copy.getDegrees();
        double up = 0;
        for(ui i = 0; i < new_vertices_count_; i++){
            up = std::max(up, (double) deg[i] / weight_[i]);
        }
        copy.subgraph_density_upper_bound = std::min(graph.subgraph_density_upper_bound, up);
        graph = copy;
    }
    
    
}

void Graph::createVertexWeightedGraph(Graph &vw_graph, double ratio) {
    Graph tmp(false, 2 * vertices_count_, true);
    double sqrt_ratio = sqrt(ratio);
    for (VertexID u = 0; u < vertices_count_; u++) {
        for (auto v: adj_[0][u]) {
            tmp.addUndirectedEdge(u, v + vertices_count_);
        }
    }
    for (VertexID u = 0; u < vertices_count_; u++) {
        tmp.weight_[u] = 1.0 / 2 / sqrt_ratio;
        tmp.weight_[u + vertices_count_] = sqrt_ratio / 2;
    }
    vw_graph = tmp;
}

void Graph::addUndirectedEdge(VertexID begin, VertexID end) {
    if (begin > adj_[0].size() || end > adj_[0].size()) {
        adj_[0].resize(std::max(begin, end) + 1);
        deg_[0].resize(std::max(begin, end) + 1, 0);
    }
    adj_[0][begin].push_back(end);
    adj_[0][end].push_back(begin);
    deg_[0][begin] += 1;
    deg_[0][end] += 1;
    edges_count_++;
};

void Graph::addDirectedEdge(VertexID begin, VertexID end) {
//    if (begin > adj_[0].size() || end > adj_[0].size()){
//        int size = std::max(int(begin), int(end)) + 1;
//        adj_[0].resize(size);
//        adj_[1].resize(size);
//        deg_[0].resize(size, 0);
//        deg_[1].resize(size, 0);
//    }
    adj_[0][begin].push_back(end);
    adj_[1][end].push_back(begin);
    deg_[0][begin] += 1;
    deg_[1][end] += 1;
    edges_count_++;
};

//void Graph::deleteEdge(VertexID begin, VertexID end) {
//    if (is_directed_) {
//        auto out_edges = adj_[0][begin];
//        out_edges.erase(std::remove(out_edges.begin(), out_edges.end(), end), out_edges.end());
//
//
//        auto &in_edges = adj_[1][end];
//        in_edges.erase(std::remove(in_edges.begin(), in_edges.end(), begin), in_edges.end());
//
//        deg_[0][begin]--;
//        deg_[1][end]--;
//
//        edges_count_--;
//    } else {
//        auto &edges_begin = adj_[0][begin];
//        edges_begin.erase(std::remove(edges_begin.begin(), edges_begin.end(), end), edges_begin.end());
//
//
//        auto &edges_end = adj_[0][end];
//        edges_end.erase(std::remove(edges_end.begin(), edges_end.end(), begin), edges_end.end());
//
//        deg_[0][begin]--;
//        deg_[0][end]--;
//
//        edges_count_--;
//    }
//};
//
//void Graph::addVertex(VertexID vertex_id) {
//    if (std::find(vertices.begin(), vertices.end(), vertex_id) != vertices.end()) {
//        std::cout << "Vertex " << vertex_id << " already exists." << std::endl;
//        return;
//    }
//
//    vertices.push_back(vertex_id);
//
//    vertices_count_++;
//
//    if (is_directed_) {
//        adj_[0].resize(vertices_count_);
//        adj_[1].resize(vertices_count_);
//        deg_[0].resize(vertices_count_, 0);
//        deg_[1].resize(vertices_count_, 0);
//    } else {
//        adj_[0].resize(vertices_count_);
//        deg_[0].resize(vertices_count_, 0);
//    }
//};

ui Graph::getEdgesCount() const {
    return edges_count_;
}

ui Graph::getVerticesCount() const {
    return vertices_count_;
}

double Graph::getMaxdeg() {
    double res = 0;
    for (ui i = 0; i < vertices_count_; i++) {
        res = std::max(res, deg_[0][i] / weight_[i]);
    }
    return res;
}

void Graph::setVerticesCount(ui n) {
    vertices_count_ = n;
    if (is_directed_) {
        for (int i = 0; i < 2; i++) {
            adj_[i].clear();
            deg_[i].clear();
            adj_[i].resize(static_cast<unsigned long>(n));
            deg_[i].resize(static_cast<unsigned long>(n));
        }
    } else {
        adj_[0].resize(static_cast<unsigned long>(n));
        weight_.resize(static_cast<double>(n));
        deg_[0].resize(static_cast<unsigned long>(n));
        vertices[0].resize(static_cast<unsigned long>(n));
    }
}

std::vector<VertexID> *Graph::getVertices() {
    return vertices;
}

std::vector<VertexID> &Graph::getNeighbors(VertexID i) const {
    return adj_[0][i];
}

std::vector<VertexID> &Graph::getOutNeighbors(VertexID i) const {
    return adj_[0][i];
}

std::vector<VertexID> &Graph::getInNeighbors(VertexID i) const {
    return adj_[1][i];
}

std::vector<std::vector<std::vector<VertexID>>> Graph::getAdjList(){
    std::vector<std::vector<std::vector<VertexID>>> adj(2);
    adj[0] = adj_[0];
    adj[1] = adj_[1];
    return adj;
}

std::vector<ui> &Graph::getDegrees() const {
    if (is_directed_) {
        printf("getDegrees() is only for undirected graphs.");
        exit(1);
    }
    return deg_[0];
}

std::vector<ui> &Graph::getInDegrees() const {
    if (!is_directed_) {
        printf("getInDegrees() is only for directed graphs.");
        exit(1);
    }
    return deg_[1];
}

std::vector<ui> &Graph::getOutDegrees() const {
    if (!is_directed_) {
        printf("getOutDegrees() is only for directed graphs.");
        exit(1);
    }
    return deg_[0];
}

