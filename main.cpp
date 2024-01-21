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
#include "app.h"

int main(int argc, char **argv) {
    Args args = Args();
    args.argsParse(argc, argv);
    //std::string graph_path = args.getOption("-path");
    std::string graph_type = args.getOption("-d");
    std::string accuracy = args.getOption("-a");
    //std::string epsilon = args.getOption("-eps");
    std::string rec_type = args.getOption("-rec"); // bool
    std::string alloc_type = args.getOption("-alloc");
    std::string ext_type = args.getOption("-ext");
    std::string ver_type = args.getOption("-ver");
    //合法性检验
    //todo
    Graph graph = Graph(graph_type != "u");
//    if () {
//        graph = Graph(false);
//    } else {
//        graph = Graph(true);
//    }
    graph.loadGraphFromFile(args.getOption("-g"));
    Reduction rec;
    Allocation alloc;
    Extraction ext;
    Verification ver;
    if (accuracy == "a") {
        //todo
    }
    else {
        Reduction rec;
        Allocation alloc;
        Extraction ext;
        Verification ver;
        if(graph_type == "u"){
            double l, r;
            ui T = 1;
            FlowNetwork flow;
            LinearProgramming lp = LinearProgramming(0, 1);
            lp.Init(graph);
            CoreApp ca = CoreApp();
            ca.Init(graph);
            l = graph.subgraph_density;
            r = graph.subgraph_density_upper_bound;
            auto vertices = new std::vector<VertexID>[1];
            bool flag = true;
            while(flag){
                T <<= 1;
                if(rec_type == "flow-exact"){
                    //todo
                }
                if(rec_type == "core-exact")
                    rec.kCoreReduction(graph, l, r);
                if(rec_type == "lp-exact"){
                    //todo
                }
                if(alloc_type == "core-app")
                    alloc.UndirectedCoreAppAllocation(graph, ca);
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
            /*
            //The generation of Ratio set needs to be refined.
            //How to combine divide-and-conquer strategy with our current framework
            //needs to be considered.
            for(int i = 1; i <= graph.getVerticesCount(); i++)
                for(int j = 1; j <= graph.getVerticesCount(); j++){
                    double ratio = i / j;
                    bool flag = true;
                    double l, r;
                    FlowNetwork flow;
                    l = graph.subgraph_density;
                    r = graph.subgraph_density_upper_bound;
                    auto vertices = new std::vector<VertexID>[2];
                    while(flag){
                        if (rec_type == "xy-core")
                            rec.xyCoreReduction(graph, ratio, l, r);

                        if (alloc_type == "flow-exact")
                            alloc.flowExactAllocation(graph, flow, ratio, l, r);
                        else if (alloc_type == "cp-exact")
                            ;

                        if (ext_type == "flow-exact")
                            ext.flowExactExtraction(graph, flow, l, r, vertices);
                        else if (ext_type == "cp-exact")
                            ;

                        if (ver_type == "flow-exact")
                            flag = ver.flowExactVerification(graph, l, r);
                        else if(ver_type == "cp-exact")
                            ;
                    }
                }
            */
        }
    }
        return 0;
};