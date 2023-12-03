//
// Created by yy on 11/30/23.
//

#include "graph.h"
#include "reduction.h"
#include "flownetwork.h"
#include "args.h"
#include "allocation.h"
#include "extraction.h"
#include "verification.h"
#include <iostream>
#include <string>
#include <fstream>

int main(int argc, char **argv) {
//    Args args;
//    args.argsParse(argc, argv);
//    std::cout << args.getOption("-a") << std::endl;
    Graph graph(true);
//    graph.loadGraphFromFile("/home/yy/DensestSubgraph/data/xycores.txt");
    graph.loadGraphFromFile("/home/yy/DensestSubgraph/data/counter_example_for_ksapp.txt");
    Allocation alloc;
    Extraction ext;
    Verification ver;
    FlowNetwork flow;
    //ratio 集合这里暂时写的比较简单
    for(int i = 1; i <= graph.getVerticesCount(); i++)
        for(int j = 1; j <= graph.getVerticesCount(); j++) {
            double ratio = i / j;
            bool flag = true;
            double l, r;
            l = graph.subgraph_density;
            r = graph.subgraph_density_upper_bound;
            auto vertices = new std::vector<VertexID>[2];

            //内层循环
            while(flag) {
                alloc.flowExactAllocation(graph, flow, ratio, l, r);
                ext.flowExactExtraction(graph, flow, l, r, vertices);
                flag = ver.flowExactVerification(graph, l, r);
            }
            printf("density %f, S/T %d/%d\n", graph.subgraph_density, graph.vertices[0].size(), graph.vertices[1].size());
        }

//    Reduction rec;
//    auto core = rec.xyCoreDecomposition(graph, 2, 2);
//    for (int i = 0; i < 2; i++) {
//        for (auto &vertex: core[i]) {
//            std::cout << vertex << " ";
//        }
//        std::cout << std::endl;
//    }
//    std::ifstream file("/home/yy/DensestSubgraph/data/maxflow.txt");
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
//    fn.getMaxFlow(0, 5);
//    fn.getMinCut(0, 5, S, T);
//
//    for(auto s: S)
//        std::cout << s << " ";



    return 0;
};
