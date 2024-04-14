#ifndef DENSESTSUBGRAPH_FLOWNETWORK_H
#define DENSESTSUBGRAPH_FLOWNETWORK_H

#include <vector>
#include <queue>
#include "types.h"
#include "graph.h"

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
//    FlowNetwork();

    FlowNetwork(ui vertices_count = 0);

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

class Edge {
public:
    ui from, to, rev;
    double flow;
    Edge(ui from, ui to, double flow, ui rev);
};
class Dinic{
public:
    std::vector<std::vector<Edge> > e;
    std::vector<int> dis;
    ui s,t;
    void Init(ui n);
    void addEdge(ui from, ui to, double flow);
    bool bfs();
    double dfs(ui u, double flow);
    double solve(double iter);
};
#endif //DENSESTSUBGRAPH_FLOWNETWORK_H

