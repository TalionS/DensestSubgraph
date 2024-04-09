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
    std::vector<ui> *deg_;

public:
   bool is_mapped;
   std::vector<VertexID> map;
//    std::vector<std::vector<VertexID>> vertex_ids;
    std::vector<double> weight_;
    std::vector<VertexID> *vertices;
    double subgraph_density_upper_bound;
    double subgraph_density_lower_bound;
    double subgraph_density;

public:
    explicit Graph(bool is_directed, ui n = 0, bool is_weighted = false);
    Graph& operator=(const Graph& other);
    ~Graph();

public:
    void init();
    void loadGraphFromFile(const std::string &dir, bool is_sample = false, double sample_rate = 1);
    void createVertexWeightedGraph(Graph &vw_graph, double ratio);
    void removeMultiEdges(Graph &graph);
    void coreOrder(Graph &graph);
    std::vector<ui> CoreOrder(Graph &graph);
    void coreReduce(Graph &graph, ui k = 0, bool weighted = false);
    void addDirectedEdge(VertexID begin, VertexID end);
    void addUndirectedEdge(VertexID begin, VertexID end);
    void sample(double p);
//    void deleteEdge(VertexID begin, VertexID end);
//    void addVertex(VertexID vertex_id);
//    void deleteVertex(VertexID vertex_id);

public:
    ui getEdgesCount() const;
    ui getVerticesCount() const;
    double getMaxdeg();
    void setVerticesCount(ui vertices_count);
    std::vector<VertexID>* getVertices();
    std::vector<VertexID> &getNeighbors(VertexID i) const;
    std::vector<VertexID> &getOutNeighbors(VertexID i) const;
    std::vector<VertexID> &getInNeighbors(VertexID i) const;
    std::vector<std::vector<std::vector<VertexID>>> getAdjList();
    std::vector<ui> &getDegrees() const;
    std::vector<ui> &getInDegrees() const;
    std::vector<ui> &getOutDegrees() const;
};


#endif //DENSESTSUBGRAPH_GRAPH_H
