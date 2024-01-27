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
    bool is_seq = args.getOption("-seq") == "t";

    double epsilon = std::stod(args.getOption("-eps"));
    double learning_rate = std::stod(args.getOption("-lr"));

    ui order_type = std::stoi(args.getOption("-o"));
    ui iter_num = std::stoi(args.getOption("-it"));
    std::string red_type = args.getOption("-red");
    std::string alloc_type = args.getOption("-alloc");
    std::string ext_type = args.getOption("-ext");
    std::string ver_type = args.getOption("-ver");
    //todo
    Graph graph = Graph(is_directed);
    graph.loadGraphFromFile(args.getOption("-path"));
    clock_t begin = clock();
    Reduction red;
    Allocation alloc;
    Extraction ext;
    Verification ver;
    if (!is_exact) {
        if (!is_directed) {
            //todo
        } else {
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
                ratio_count++;
                bool flag = true;
                bool is_init_red = false;
                bool is_init_lp = false;
                bool is_reduced = false;
                bool is_stable_set = false;
                double rho, rho_c;
                double l = learning_rate * graph.subgraph_density;
                double r = graph.subgraph_density_upper_bound;
                FlowNetwork flow;
                LinearProgramming lp(is_directed, 0, 0, 0, order_type);
                std::vector<std::pair<VertexID, VertexID>> edges;
                std::vector<std::vector<VertexID>> vertices(2);
                std::pair<ui, ui> best_pos(0, 0);
                Graph subgraph(is_directed, graph.getVerticesCount());

                while (flag) {
                    if (!is_reduced || alloc_type != "CP") {
                        is_reduced = true;
                        if (red_type == "exact-xy-core") {
                            red.xyCoreReduction(graph, subgraph, ratio, l, r, is_init_red,
                                                is_dc, false, true);
                        } else if (red_type == "appro-xy-core") {
                            red.xyCoreReduction(graph, subgraph, ratio, l, r, is_init_red,
                                                is_dc, false, false);
                        }
                    }
                    if (alloc_type == "CP")
                        alloc.directedCPAllocation(subgraph, lp, iter_num, is_init_lp, ratio, !is_seq, is_exp);

                    if (ext_type == "CP")
                        ext.directedCPExtraction(subgraph, lp, best_pos, vertices, ratio, ratio_o, ratio_p, rho, rho_c);

                    if (ver_type == "CP")
                        flag = ver.directedCPVerification(graph, subgraph, lp, best_pos, vertices, ratio, rho, rho_c,
                                                          ratio_o, ratio_p, is_stable_set, edges, epsilon);
                }
            }
            printf("ratio count: %d, density: %f, S/T: %d/%d\n", ratio_count, graph.subgraph_density,
                   graph.vertices[0].size(), graph.vertices[1].size());
        }
    } else {
//        Reduction red;
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
                    red.kCoreReduction(graph, l, r);
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
                                                  0)) {
                ratio_count++;
                bool flag = true;
                bool is_init_red = false;
                bool is_init_lp = false;
                bool is_reduced = false;
                bool is_stable_set = false;
                double rho, rho_c;
                double l = learning_rate * graph.subgraph_density;
                double r = graph.subgraph_density_upper_bound;
                FlowNetwork flow;
                LinearProgramming lp(is_directed, 0, 0, 0, order_type);
                std::vector<std::pair<VertexID, VertexID>> edges;
                std::vector<std::vector<VertexID>> vertices(2);
                std::pair<ui, ui> best_pos(0, 0);
                Graph subgraph(is_directed, graph.getVerticesCount());

                while (flag) {
                    if (!is_reduced || alloc_type != "CP") {
                        is_reduced = true;
                        if (red_type == "exact-xy-core") {
                            red.xyCoreReduction(graph, subgraph, ratio, l, r, is_init_red,
                                                is_dc, false, true);
                        } else if (red_type == "appro-xy-core") {
                            red.xyCoreReduction(graph, subgraph, ratio, l, r, is_init_red,
                                                is_dc, false, false);
                        }
                    }
                    if (alloc_type == "CP")
                        red.stableSetReduction(subgraph, lp, edges, is_stable_set);

                    if (alloc_type == "CP")
                        alloc.directedCPAllocation(subgraph, lp, iter_num, is_init_lp, ratio, !is_seq, is_exp);
                    if (alloc_type == "flow-exact")
                        alloc.flowExactAllocation(subgraph, flow, ratio, l, r, is_dc);

                    if (ext_type == "CP")
                        ext.directedCPExtraction(subgraph, lp, best_pos, vertices, ratio, ratio_o, ratio_p, rho, rho_c);
                    if (ext_type == "flow-exact")
                        ext.flowExactExtraction(graph, ratio, flow, l, r, ratio_o, ratio_p);

                    if (ver_type == "CP")
                        flag = ver.directedCPVerification(graph, subgraph, lp, best_pos, vertices, ratio, rho, rho_c,
                                                          ratio_o, ratio_p, is_stable_set, edges, 0);
                    if (ver_type == "flow-exact")
                        flag = ver.flowExactVerification(graph, l, r);
                }
            }
            printf("ratio count: %d, density: %f, S/T: %d/%d\n", ++ratio_count, graph.subgraph_density,
                   graph.vertices[0].size(), graph.vertices[1].size());
        }

    }
    clock_t end = clock();
    printf("time: %f", (double) (end - begin) / CLOCKS_PER_SEC);
    return 0;
};