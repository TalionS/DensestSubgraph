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
#include "wcore.h"
#include <boost/heap/fibonacci_heap.hpp>
using Heap = boost::heap::fibonacci_heap<std::pair<int, VertexID>>;

int main(int argc, char **argv) {
//    Args args;
//    args.argsParse(argc, argv);
//    std::cout << args.getOption("-a") << std::endl;
    clock_t begin = clock();
    Graph graph(true);
//    graph.loadGraphFromFile("../data/counter_example_for_ksapp.txt");
//    graph.loadGraphFromFile("../data/xycores.txt");
//    graph.loadGraphFromFile("../data/MI.txt");
//    density 7.606087, S/T 13/12
    graph.loadGraphFromFile("../data/AD1.txt");
//    density 31.681085, S/T 453/195
    auto ratioSelect = RatioSelection(graph, true);
    Reduction red;
    Allocation alloc;
    Extraction ext;
    Verification ver;
    FlowNetwork flow;
    XYCore xycore;
    std::pair<ui, ui> max_core_num;
//    alloc.coreApproAllocation(graph, max_core_num);
//    xycore.xyCoreInitialization(graph);
//    Graph x_y_core = Graph(true, graph.getVerticesCount());
//    xycore.generateXYCore(graph, x_y_core, 2, 14);
//    printf("%d\n", x_y_core.getEdgesCount());
//    WCore w_core;
//    red.wCoreReduction(x_y_core, w_core);
//    w_core.getMaxCNPair(x_y_core, max_core_num);
    ui ratio_count = 0;
    std::vector<std::vector<VertexID>> vertices(2);
    bool is_dc = true;
    bool is_init_ratio = false;
    bool is_init_red = false;
    auto ratio = std::pair<double, double>(0, 0);
    double ratio_o, ratio_p;
    ui s_size = 0, t_size = 0;
    std::vector<std::pair<VertexID, VertexID>> edges;
    graph.subgraph_density = 0;
    //ratio 集合这里暂时写的比较简单
    while(ratioSelect.ratioSelection(
            graph.getVerticesCount(),
            ratio,
            is_init_ratio, false,
            false,
            false,
//            true,
//            false,
            s_size,
            t_size,
            ratio_o,
            ratio_p,
            graph.subgraph_density,
            0))
    {
        s_size = 0;
        t_size = 0;
        bool flag = true;
        bool is_init_lp = false;
        double l, r;
        l = 1 * graph.subgraph_density;
//        l = 0;
        r = sqrt(graph.subgraph_density_upper_bound);
        LinearProgramming lp = LinearProgramming(true, 1);
        ui T = 2;
        bool is_core = false;
        std::pair<ui, ui> best_pos(0, 0);
        double rho, rho_c;
        Graph x_y_core = Graph(true, graph.getVerticesCount());
        bool stable_set_reduction = false;
//        Graph x_y_core = graph;
        //内层循环
        ui edges_count = 0;
        std::vector<std::vector<ui>> degrees(2);
        std::vector<Heap> heap;
        std::vector<std::vector<Heap::handle_type>> handles;
        std::vector<std::vector<bool>> is_peeled;
        std::vector<ui> vertices_count;
        ui cur = 0;
        while(flag) {
//            Graph x_y_core = Graph(true, graph.getVerticesCount());
            if (!is_core)
                red.xyCoreReduction(graph, x_y_core, ratio, l, r, is_init_red, is_dc);
            is_core = true;
            red.stableSetReduction(x_y_core, lp, edges, stable_set_reduction);
//            T += 100;
            T <<= 1;
            printf("%d\n", T);
//            alloc.directedCPAllocation(graph, lp, T, is_init_lp, ratio, true);
//            ext.directedVWApproExtraction(graph, lp, vertices, ratio, rho, rho_c);
//            flag = ver.directedVWApproVerification(graph, lp, vertices, rho, rho_c, 0);


            alloc.directedCPAllocation(x_y_core, lp, T, is_init_lp, ratio);
            ext.directedCPExtraction(x_y_core, lp, best_pos, vertices, ratio, ratio_o, ratio_p, rho, rho_c);
            flag = ver.directedCPVerification(graph, x_y_core, lp, best_pos, vertices, ratio, rho, rho_c, ratio_o,
                                              ratio_p, stable_set_reduction,
                                              edges, 0);

//            alloc.directedPMApproAllocation(graph, ratio, 1, edges_count, vertices, degrees, is_init_lp);
//            ext.directedPMApproExtraction(graph, edges_count, vertices);
//            flag = ver.directedPMApproVerification(vertices);

//            alloc.directedFixedKSApproAllocation(graph, ratio, cur, heap, handles, is_peeled, edges_count);
//            alloc.directedBSApproAllocation(graph, ratio, heap, handles, is_peeled, edges_count, vertices_count, is_init_lp);
//            alloc.directedKSApproAllocation(graph, heap, handles, is_peeled, edges_count, is_init_lp);
//            ext.directedBSApproExtraction(graph, is_peeled, vertices);
//            flag = ver.directedFixedKSApproVerification(graph, cur, edges_count, vertices);
//            flag = ver.directedBSApproVerification(graph, edges_count, vertices);

//            alloc.flowExactAllocation(x_y_core, flow, ratio, l, r, is_dc);
//            ext.flowExactExtraction(graph, flow, l, r, s_size, t_size);
//            flag = ver.flowExactVerification(graph, l, r);
//            WCore w_core;
//            red.wCoreReduction(graph, w_core);
//            w_core.getMaxCNPair(graph, max_core_num);
//            ext
        }
        printf("ratio_count %d, ratio %f, density %f, S/T %d/%d\n", ++ratio_count, ratio.first / ratio.second, graph.subgraph_density, graph.vertices[0].size(), graph.vertices[1].size());
//        printf("ratio_count %d, sqrt_ratio %.4f, density %f, S/T %d/%d\n", ++ratio_count, sqrt(ratio.first / ratio.second), graph.subgraph_density, graph.vertices[0].size(), graph.vertices[1].size());
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
//    std::ifstream file("../data/maxflow2.txt");
//    ui n;
//    file >> n;
//    FlowNetwork fn(n);
//    VertexID from, to;
//    double cap;
//    while (file >> from >> to >> cap){
//        fn.addEdge(from, to, cap);
//    }
//    file.close();
//
//    std::vector<VertexID> S, T;
//    auto maxflow = fn.getMaxFlow(0, n - 1);
//    printf("%f", maxflow);
//    fn.getMinCut(0, 5, S, T);
//
//    for(auto s: S)
//        std::cout << s << " ";



    return 0;
};
