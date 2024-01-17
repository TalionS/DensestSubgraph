//
// Created by yy on 23-11-13.
//

#include "graph.h"
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
            vertices[i].resize(static_cast<unsigned long>(n));
        }
    } else {
        adj_ = new std::vector<std::vector<VertexID>>[1];
        deg_ = new std::vector<ui>[1];
        vertices = new std::vector<VertexID>[1];
        adj_[0].resize(static_cast<unsigned long>(n));
        weight_.resize(static_cast<unsigned long>(n));
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

        // Resize and copy adjacency lists, degree vectors, and vertices vectors
        if (is_directed_) {
            adj_[0] = other.adj_[0];
            adj_[1] = other.adj_[1];
            deg_[0] = other.deg_[0];
            deg_[1] = other.deg_[1];
            vertices[0] = other.vertices[0];
            vertices[1] = other.vertices[1];
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

void Graph::loadGraphFromFile(const std::string &dir) {
    std::ifstream file(dir);
    if (!file.is_open()) {
        std::cout << "Unable to open file " << dir << "." << std::endl;
        exit(-1);
    }

    ui edges_count;
    file >> vertices_count_ >> edges_count;
    if (is_directed_) {
        for (int i = 0; i < 2; i++) {
            adj_[i].resize(static_cast<VertexID>(vertices_count_));
            deg_[i].resize(static_cast<VertexID>(vertices_count_), 0);
        }
        VertexID begin, end;
        while (file >> begin >> end) {
            addDirectedEdge(begin, end);
        }
    } else {
        adj_[0].resize(static_cast<VertexID>(vertices_count_));
        deg_[0].resize(static_cast<VertexID>(vertices_count_), 0);
        VertexID begin, end;
        while (file >> begin >> end) {
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
};

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

ui Graph::getMaxdeg() {
    ui res = 0;
    for (ui i = 0; i < vertices_count_; i++) {
        res = std::max(res, (ui) deg_[i].size());
    }
    return res;
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

