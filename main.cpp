//
// Created by yy on 23-11-14.
//
#include "graph.h"
#include "reduction.h"
#include "flownetwork.h"
#include <iostream>
#include <fstream>

int main() {
    Graph obj(true);
    obj.loadGraphFromFile("/home/yy/DensestSubgraph/data/xycores.txt");
    Reduction rec;
    auto core = rec.xyCoreDecomposition(obj, 2, 2);
    for (int i = 0; i < 2; i++) {
        for(auto &vertex: core[i]){
            std::cout << vertex << " ";
        }
        std::cout << std::endl;
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