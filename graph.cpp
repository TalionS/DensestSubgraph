//
// Created by yy on 23-11-13.
//

#include "graph.h"
#include <fstream>
#include <iostream>
#include <assert.h>

Graph::Graph(bool is_directed) : is_directed_(is_directed), vertices_count_(0), edges_count_(0){
    if (is_directed_) {
        adj_ = new std::vector<std::vector<VertexID>>[2];
        deg_ = new std::vector<ui>[2];
    } else {
        adj_ = new std::vector<std::vector<VertexID>>[1];
        deg_ = new std::vector<ui>[1];
    }
}

Graph::~Graph() {
    delete[] adj_;
    delete[] deg_;
}

void Graph::loadGraphFromFile(const char *dir){
    std::ifstream file(dir);
    if(!file.is_open()){
        std::cout << "Unable to open file "<< dir<< "."<< std::endl;
        exit(-1);
    }

    ui edges_count;
    file >> vertices_count_ >> edges_count;
    vertices_.resize(static_cast<VertexID>(vertices_count_), 0);
    for(VertexID i = 0; i < static_cast<VertexID>(vertices_count_); i++){
        vertices_[i] = i;
    }
    if(is_directed_) {
        for(int i = 0; i < 2; i++){
            adj_[i].resize(static_cast<VertexID>(vertices_count_));
            deg_[i].resize(static_cast<VertexID>(vertices_count_), 0);
        }
        VertexID begin, end;
        while(file >> begin >> end){
            addDirectedEdge(begin, end);
        }
    }
    else{
        adj_[0].resize(static_cast<VertexID>(vertices_count_));
        deg_[0].resize(static_cast<VertexID>(vertices_count_), 0);
        VertexID begin, end;
        while(file >> begin >> end){
            addUndirectedEdge(begin, end);
        }
    }

    assert(edges_count_==edges_count);
    file.close();
};

void Graph::addUndirectedEdge(VertexID begin, VertexID end){
    adj_[0][begin].push_back(end);
    adj_[0][end].push_back(begin);
    deg_[0][begin] += 1;
    deg_[0][end] += 1;
    edges_count_++;
};

void Graph::addDirectedEdge(VertexID begin, VertexID end){
    adj_[0][begin].push_back(end);
    adj_[1][end].push_back(begin);
    deg_[0][begin] += 1;
    deg_[1][end] += 1;
    edges_count_++;
};

void Graph::deleteEdge(VertexID begin, VertexID end){
    if(is_directed_) {
        auto out_edges = adj_[0][begin];
        out_edges.erase(std::remove(out_edges.begin(), out_edges.end(), end), out_edges.end());


        auto& in_edges = adj_[1][end];
        in_edges.erase(std::remove(in_edges.begin(), in_edges.end(), begin), in_edges.end());

        deg_[0][begin]--;
        deg_[1][end]--;

        edges_count_--;
    }
    else {
        auto& edges_begin = adj_[0][begin];
        edges_begin.erase(std::remove(edges_begin.begin(), edges_begin.end(), end), edges_begin.end());


        auto& edges_end = adj_[0][end];
        edges_end.erase(std::remove(edges_end.begin(), edges_end.end(), begin), edges_end.end());

        deg_[0][begin]--;
        deg_[0][end]--;

        edges_count_--;
    }
};

void Graph::addVertex(VertexID vertex_id){
    if (std::find(vertices_.begin(), vertices_.end(), vertex_id) != vertices_.end()) {
        std::cout << "Vertex " << vertex_id << " already exists." << std::endl;
        return;
    }

    vertices_.push_back(vertex_id);

    vertices_count_++;

    if(is_directed_) {
        adj_[0].resize(vertices_count_);
        adj_[1].resize(vertices_count_);
        deg_[0].resize(vertices_count_, 0);
        deg_[1].resize(vertices_count_, 0);
    } else {
        adj_[0].resize(vertices_count_);
        deg_[0].resize(vertices_count_, 0);
    }
};

const ui Graph::getEdgesCount() const {
    return edges_count_;
}

const ui Graph::getVerticesCount() const {
    return vertices_count_;
}

const std::vector<VertexID>& Graph::getNeighbors(VertexID i) const {
    return adj_[0][i];
}

const std::vector<VertexID>& Graph::getOutNeighbors(VertexID i) const {
    return adj_[0][i];
}

const std::vector<VertexID>& Graph::getInNeighbors(VertexID i) const {
    return adj_[1][i];
}

const std::vector<ui>& Graph::getDegrees() const {
    if(is_directed_) {
        printf("getDegrees() is only for undirected graphs.");
        exit(1);
    }
    return deg_[0];
}

const std::vector<ui>& Graph::getInDegrees() const {
    if(!is_directed_) {
        printf("getInDegrees() is only for directed graphs.");
        exit(1);
    }
    return deg_[1];
}

const std::vector<ui>& Graph::getOutDegrees() const {
    if(!is_directed_) {
        printf("getOutDegrees() is only for directed graphs.");
        exit(1);
    }
    return deg_[0];
}

