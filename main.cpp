//
// Created by yy on 23-11-14.
//
#include "graph.h"
#include "reduction.h"
#include "flownetwork.h"
#include <iostream>
#include <fstream>

int main(){
    Graph obj(false);
    obj.loadGraphFromFile("/home/yy/DensestSubgraph/data/cores.txt");
    std::vector<VertexID> a = obj.getNeighbors(static_cast<VertexID>(0));
    Reduction rec;
    std::vector<ui> core;
    core = rec.coreDecomposition(obj);
    for(int i=0; i < obj.getVerticesCount(); i++){
        std::cout << core[i] << " ";
    }
//    std::ifstream file("/home/yy/DensestSubgraph/data/maxflow2.txt");
//    ui n;
//    file >> n;
//    FlowNetwork fn(n);
//    VertexID from, to;
//    ui cap;
//    while (file >> from >> to >> cap){
//        fn.addEdge(from, to, cap);
//    }
//    file.close();
//
//    std::vector<VertexID> S, T;
//    std::cout << fn.getMinCut(0, 3, S, T)<< std::endl;
//
//    for(auto s: S)
//        std::cout << s << " ";



    return 0;
};