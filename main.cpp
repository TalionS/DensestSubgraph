//
// Created by yy on 23-11-14.
//
#include "graph.h"
#include "reduction.h"
#include "allocation.h"
#include "extraction.h"
#include "verification.h"
#include "flownetwork.h"
#include "lp.h"
#include <iostream>
#include <fstream>
#include "args.h"
#include "ratioselection.h"
#include "app.h"
#include <ctime>

int main(int argc, char **argv) {
    Args args = Args();
    args.argsParse(argc, argv);
    bool is_directed = args.getOption("-t") == "d";
    bool is_exact = args.getOption("-a") == "e";
    bool is_vw = args.getOption("-vw") == "t";
    bool is_parallel = args.getOption("-p") == "t";
    bool is_exp = args.getOption("-exp") == "t";
    bool is_dc = args.getOption("-dc") == "t";
    double epsilon = std::stod(args.getOption("-eps"));
    double learning_rate = std::stod(args.getOption("-lr"));
    ui iter_num = std::stoi(args.getOption("-it"));
    std::string red_type = args.getOption("-red");
    std::string alloc_type = args.getOption("-alloc");
    std::string ext_type = args.getOption("-ext");
    std::string ver_type = args.getOption("-ver");
    std::string order_type = args.getOption("-o");
    std::string update_type = args.getOption("-s");
    //todo
    Graph graph = Graph(is_directed);
    graph.loadGraphFromFile(args.getOption("-path"));
    clock_t begin = clock();
    Reduction rec;
    Allocation alloc;
    Extraction ext;
    Verification ver;
    if (!is_exact) {
        //todo
    } else {
//        Reduction rec;
//        Allocation alloc;
//        Extraction ext;
//        Verification ver;
        if (!is_directed) {
            double l, r;
            ui T = 1;
            FlowNetwork flow;
            LinearProgramming lp = LinearProgramming(0, 1);
            lp.Init(graph);
//            CoreApp ca = CoreApp();
//            ca.Init(graph);
            l = graph.subgraph_density;
            r = graph.subgraph_density_upper_bound;
            auto vertices = new std::vector<VertexID>[1];
            bool flag = true;
            while (flag) {
                T <<= 1;
                if (red_type == "flow-exact") {
                    //todo
                }
                if (red_type == "core-exact")
                    rec.kCoreReduction(graph, l, r);
                if (red_type == "lp-exact") {
                    //todo
                }
//                if (alloc_type == "core-app")
//                    alloc.UndirectedCoreAppAllocation(graph, ca);
                if (alloc_type == "flow-exact")
                    alloc.UndirectedflowExactAllocation(graph, flow, l, r);
                if (alloc_type == "lp-exact")
                    alloc.UndirectedlpAllocation(graph, lp, T);
                if (ext_type == "flow-exact")
                    ext.UndirectedflowExactExtraction(graph, flow, l, r, vertices);
                if (ext_type == "lp-exact")
                    ext.UndirectedlpExactExtraction(graph, lp, vertices);
                if (ver_type == "flow-exact")
                    flag = ver.UndirectedflowExactVerification(graph, l, r);
                if (ver_type == "lp-exact")
                    flag = ver.UndirectedlpVerification(graph, lp, flow, vertices);
            }
            //todo
        } else {
            //The generation of Ratio set needs to be refined.
            //How to combine divide-and-conquer strategy with our current framework
            //needs to be considered.
            std::pair<double, double> ratio;
            double ratio_o, ratio_p;
            RatioSelection ratio_selection(graph);
            bool is_init_ratio = false;
            ui ratio_count = 0;
            while (ratio_selection.ratioSelection(graph.getVerticesCount(),
                                                  ratio,
                                                  is_init_ratio,
                                                  is_vw,
                                                  is_dc,
                                                  ratio_o,
                                                  ratio_p,
                                                  graph.subgraph_density,
                                                  epsilon)) {

                bool flag = true;
                bool is_init_red = false;
                double l, r;
                FlowNetwork flow;
                l = learning_rate * graph.subgraph_density;
                r = graph.subgraph_density_upper_bound;
                auto vertices = new std::vector<VertexID>[2];
                while (flag) {
                    Graph x_y_core(is_directed, graph.getVerticesCount());
                    if (red_type == "exact-xy-core") {
                        rec.xyCoreReduction(graph, x_y_core, ratio, l, r, is_init_red,
                                            is_dc, false, true);
                    } else if (red_type == "appro-xy-core") {
                        rec.xyCoreReduction(graph, x_y_core, ratio, l, r, is_init_red,
                                            is_dc, false, false);
                    }

                    if (alloc_type == "flow-exact")
                        alloc.flowExactAllocation(graph, flow, ratio, l, r, is_dc);

                    if (ext_type == "flow-exact")
                        ext.flowExactExtraction(graph, ratio, flow, l, r, ratio_o, ratio_p);

                    if (ver_type == "flow-exact")
                        flag = ver.flowExactVerification(graph, l, r, <#initializer#>, 0, <#initializer#>);
                }
                printf("ratio count: %d, density: %f, S/T: %d/%d\n", ++ratio_count, graph.subgraph_density, graph.vertices[0].size(), graph.vertices[1].size());
            }
        }
    }
    clock_t end = clock();
    printf("time: %f", (double) (end - begin) / CLOCKS_PER_SEC);
    return 0;
};