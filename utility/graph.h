//
// Created by yy on 23-11-13.
//

#ifndef DENSESTSUBGRAPH_GRAPH_H
#define DENSESTSUBGRAPH_GRAPH_H

#include <vector>
#include <iostream>
#include "types.h"
#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>

class Graph {
private:
    bool is_directed_;
    bool is_weighted_;
    ui vertices_count_;
    ui edges_count_;

    std::vector<std::vector<VertexID>> *adj_;
    std::vector<std::vector<ui>> weight_;
    std::vector<ui> *deg_;

public:
    std::vector<VertexID> *vertices;
    double subgraph_density_upper_bound;
    double subgraph_density_lower_bound;
    double subgraph_density;

public:
    explicit Graph(bool is_directed, ui n = 0, bool is_weighted = false);
    Graph& operator=(const Graph& other);
    ~Graph();

public:
    void loadGraphFromFile(const std::string &dir);
    void addDirectedEdge(VertexID begin, VertexID end);
    void addUndirectedEdge(VertexID begin, VertexID end, ui weight = 1);

//    void deleteEdge(VertexID begin, VertexID end);
//    void addVertex(VertexID vertex_id);
//    void deleteVertex(VertexID vertex_id);

public:
    ui getEdgesCount() const;
    ui getVerticesCount() const;
    ui getMaxdeg();
    std::vector<VertexID>* getVertices();
    std::vector<VertexID> &getNeighbors(VertexID i) const;
    std::vector<VertexID> &getOutNeighbors(VertexID i) const;
    std::vector<VertexID> &getInNeighbors(VertexID i) const;
    std::vector<ui> &getDegrees() const;
    std::vector<std::vector<ui>> getWeights() const;
    std::vector<ui> &getInDegrees() const;
    std::vector<ui> &getOutDegrees() const;
};


#endif //DENSESTSUBGRAPH_GRAPH_H
