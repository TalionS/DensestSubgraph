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

int main(int argc, char **argv) {
    Args args = Args();
    args.argsParse(argc, argv);
    bool is_directed =  args.getOption("-t") == "d";
    bool is_exact = args.getOption("-a") == "e";
    bool is_vw = args.getOption("-vw") == "t";
    bool is_parallel = args.getOption("-p") == "t";
    double epsilon = std::stod(args.getOption("-eps"));
    double learning_rate = std::stod(args.getOption("-lr"));
    std::string red_type = args.getOption("-red");
    std::string alloc_type = args.getOption("-alloc");
    std::string ext_type = args.getOption("-ext");
    std::string ver_type = args.getOption("-ver");
    std::string order_type = args.getOption("-o");
    std::string update_type = args.getOption("-s");
    //todo
    Graph graph = Graph(args.getOption("-t") != "u");
    graph.loadGraphFromFile(args.getOption("-path"));
    Reduction rec;
    Allocation alloc;
    Extraction ext;
    Verification ver;
    if (!is_exact) {
        //todo
    }
    else {
//        Reduction rec;
//        Allocation alloc;
//        Extraction ext;
//        Verification ver;
        if(!is_directed){
            double l, r;
            ui T = 1;
            FlowNetwork flow;
            LinearProgramming lp = LinearProgramming(0, 1);
            lp.Init(graph);
            l = graph.subgraph_density;
            r = graph.subgraph_density_upper_bound;
            auto vertices = new std::vector<VertexID>[1];
            bool flag = true;
            while(flag){
                T <<= 1;
                if(red_type == "flow-exact"){
                    //todo
                }
                if(red_type == "core-exact")
                    rec.kCoreReduction(graph, l, r);
                if(red_type == "lp-exact"){
                    //todo
                }
                if(alloc_type == "flow-exact")
                    alloc.UndirectedflowExactAllocation(graph, flow, l, r);
                if(alloc_type == "lp-exact")
                    alloc.UndirectedlpAllocation(graph, lp, T);
                if(ext_type == "flow-exact")
                    ext.UndirectedflowExactExtraction(graph, flow, l, r, vertices);
                if(ext_type == "lp-exact")
                    ext.UndirectedlpExactExtraction(graph, lp, vertices);
                if(ver_type == "flow-exact")
                    flag = ver.UndirectedflowExactVerification(graph, l, r);
                if(ver_type == "lp-exact")
                    flag = ver.UndirectedlpVerification(graph, lp, flow, vertices);
            }
            //todo
        }
        else{
            //The generation of Ratio set needs to be refined.
            //How to combine divide-and-conquer strategy with our current framework
            //needs to be considered.
//            std::pair<double, double> ratio;
//            RatioSelection ratio_selection(graph, true);
//            bool is_init_ratio = false;
//            bool is_vw, is_dc;
//            while(ratio_selection.ratioSelection(graph.getVerticesCount(),
//                                                 ratio,
//                                                 is_init_ratio,
//                                                 is_dc,
//                                                 is_dc))
//
//                    bool flag = true;
//                    double l, r;
//                    FlowNetwork flow;
//                    l = graph.subgraph_density;
//                    r = graph.subgraph_density_upper_bound;
//                    auto vertices = new std::vector<VertexID>[2];
//                    while(flag){
//                        if (red_type == "xy-core")
//                            rec.xyCoreReduction(graph, ratio, l, r, <#initializer#>, <#initializer#>, false, false);
//
//                        if (alloc_type == "flow-exact")
//                            alloc.flowExactAllocation(graph, flow, ratio, l, r);
//                        else if (alloc_type == "cp-exact")
//                            ;
//
//                        if (ext_type == "flow-exact")
//                            ext.flowExactExtraction(graph, flow, l, r, vertices);
//                        else if (ext_type == "cp-exact")
//                            ;
//
//                        if (ver_type == "flow-exact")
//                            flag = ver.flowExactVerification(graph, l, r);
//                        else if(ver_type == "cp-exact")
//                            ;
//                    }
//                }
        }
    }
        return 0;
};