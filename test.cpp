//
// Created by yy on 11/30/23.
//

#include "utility/graph.h"
#include "reduction/reduction.h"
#include "utility/flownetwork.h"
#include "utility/args.h"
#include "allocation/allocation.h"
#include "extraction/extraction.h"
#include "verification/verification.h"
#include "utility/ratioselection.h"
#include "utility/xycore.h"
#include <iostream>
#include <string>
#include <fstream>
#include <time.h>

int main(int argc, char **argv) {
//    Args args;
//    args.argsParse(argc, argv);
//    std::cout << args.getOption("-a") << std::endl;
    clock_t begin = clock();
    Graph graph(true);
//    graph.loadGraphFromFile("/home/yy/DensestSubgraph/data/counter_example_for_ksapp.txt");
    graph.loadGraphFromFile("/home/yy/DensestSubgraph/data/MI.txt");
    RatioSelection ratioSelect;
    Reduction red;
    Allocation alloc;
    Extraction ext;
    Verification ver;
    FlowNetwork flow;
//    XYCore xycore;
//    std::pair<ui, ui> max_core_num;
//    alloc.coreApproAllocation(graph, max_core_num);
//    xycore.xyCoreInitialization(graph);
//    printf("%d\n", xycore.getDelta(graph));
    ui ratio_count = 0;
    std::vector<std::vector<VertexID>> vertices(2);
    bool is_dc = true;
    bool is_init_ratio = false;
    bool is_init_red = false;
    auto ratio = std::pair<double, double>(0, 0);
    double ratio_o, ratio_p;
    ui s_size = 0, t_size = 0;
    //ratio 集合这里暂时写的比较简单
    while(ratioSelect.ratioSelection(
            graph.getVerticesCount(),
            ratio,
            is_init_ratio,
            false,
            true,
            s_size,
            t_size,
            ratio_o,
            ratio_p))
    {
        s_size = 0;
        t_size = 0;
        bool flag = true;
        bool is_init_lp = false;
        double l, r;
        l = graph.subgraph_density;
//        l = 0;
        r = graph.subgraph_density_upper_bound;
        LinearProgramming lp = LinearProgramming(true, 0);
        ui T = 1;
        std::pair<ui, ui> best_pos;
        double rho_c;

        //内层循环
        while(flag) {

            Graph x_y_core = Graph(true, graph.getVerticesCount());
            red.xyCoreReduction(graph, x_y_core, ratio, l, r, is_init_red, is_dc);
            T <<= 1;
            printf("%d\n", T);
            alloc.directedLPExactAllocation(x_y_core, lp, T, is_init_lp, ratio);
            ext.directedLPExactExtraction(x_y_core, lp, best_pos, vertices, ratio, ratio_o, ratio_p, rho_c);
            flag = ver.directedLPExactVerification(graph, x_y_core, lp, best_pos, vertices, ratio, rho_c);
//            alloc.flowExactAllocation(x_y_core, flow, ratio, l, r, is_dc);
//            ext.flowExactExtraction(graph, flow, l, r, s_size, t_size);
//            flag = ver.flowExactVerification(graph, red.xycore, l, r);
        }
        printf("ratio_count %d, ratio (%f, %f), density %f, S/T %d/%d\n", ++ratio_count, ratio.first, ratio.second, graph.subgraph_density, graph.vertices[0].size(), graph.vertices[1].size());
    }
    clock_t end = clock();
    printf("time: %f", (double )(end - begin) / CLOCKS_PER_SEC);

//    Reduction red;
//    auto core = red.xyCoreDecomposition(graph, 2, 2);
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
