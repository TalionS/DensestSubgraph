//
// Created by yy on 23-11-14.
//
#include "graph.h"
#include "reduction.h"
#include "allocation.h"
#include "extraction.h"
#include "verification.h"
#include "flownetwork.h"
#include <iostream>
#include <fstream>
#include "args.h"

int main(int argc, char **argv) {
    Args args = Args();
    args.argsParse(argc, argv);
    std::string graph_path = args.getOption("-path");
    std::string graph_type = args.getOption("-g_type");
    std::string accuracy = args.getOption("-appro");
    std::string epsilon = args.getOption("-eps");
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
            if(rec_type == "flow-exact"){
                //todo
            }
            if(rec_type == "core-exact"){
                //todo
            }
            if(rec_type == "lp-exact"){
                //todo
            }
            //todo
            if(alloc_type == "lp-exact"){
                //todo
            }
            //todo
            if(ext_type == "lp-exact"){
                //todo
            }
            //todo
        }
        else{

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
        }
    }
        return 0;
};