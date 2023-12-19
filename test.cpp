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
#include "ratioselection.h"
#include "xycore.h"
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
    RatioSelection ratioSelect;
    Reduction rec;
    Allocation alloc;
    Extraction ext;
    Verification ver;
    FlowNetwork flow;
    XYCore xycore;
    std::pair<ui, ui> max_core_num;
    alloc.coreApproAllocation(graph, max_core_num);
    xycore.xyCoreInitialization(graph);
    printf("%d\n", xycore.getDelta(graph));
    ui ratio_count = 0;
    bool is_dc = true;
    bool is_init_ratio = false;
    bool is_init_rec = false;
    auto vertices = new std::vector<VertexID>[2];
    auto ratio = std::pair<double, double>(0, 0);
    //ratio 集合这里暂时写的比较简单
    while(ratioSelect.ratioSelection(
            graph.getVerticesCount(),
            ratio,
            is_init_ratio,
            is_dc,
            false,
            vertices[0].size(),
            vertices[1].size()))
    {
        bool flag = true;
        double l, r;
        l = graph.subgraph_density;
        r = graph.subgraph_density_upper_bound;
        vertices = new std::vector<VertexID>[2];

        //内层循环
        while(flag) {
            Graph x_y_core = Graph(true, graph.getVerticesCount());
            rec.xyCoreReduction(graph, x_y_core, ratio, l, r, is_init_rec, is_dc);
            alloc.flowExactAllocation(graph, x_y_core, flow, ratio, l, r, is_dc);
            ext.flowExactExtraction(graph, flow, l, r, vertices);
            flag = ver.flowExactVerification(graph, l, r);
        }
        printf("ratio_count %d, density %f, S/T %d/%d\n", ++ratio_count, graph.subgraph_density, graph.vertices[0].size(), graph.vertices[1].size());
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
