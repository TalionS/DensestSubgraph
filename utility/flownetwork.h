//
// Created by yy on 11/19/23.
//

#ifndef DENSESTSUBGRAPH_FLOWNETWORK_H
#define DENSESTSUBGRAPH_FLOWNETWORK_H

#include <vector>
#include <queue>
#include "types.h"

class FlowEdge {
public:
    VertexID from, to, index;
    double capacity, flow;

    FlowEdge(VertexID from, VertexID to, double capacity, double flow, VertexID index);
};


class FlowNetwork {
public:
    std::vector<std::vector<FlowEdge>> adj_;
    std::vector<double> excess_;
    std::vector<ui> dist_;
    std::vector<ui> count_;
    std::vector<bool> active_;
    std::vector<std::vector<ui>> height_bucket_;
    std::vector<ui> nums_;
    std::vector<VertexID> ori_id_;
    std::queue<int> node_queue_;
    std::vector<VertexID> vertices_;
    ui nodes_count_;
    ui highest_active_height_;


public:
    FlowNetwork();

    FlowNetwork(ui vertices_count);

    void addEdge(VertexID from, VertexID to, double capacity);

    double getMaxFlow(VertexID s, VertexID t);

    void getMinCut(VertexID s, VertexID t, std::vector<VertexID> &S);
    

private:
    void enqueue(VertexID v);

    void push(FlowEdge &e);

    void gap(VertexID k);

    void relabel(VertexID v);

    void discharge(VertexID v);
};

#endif //DENSESTSUBGRAPH_FLOWNETWORK_H
