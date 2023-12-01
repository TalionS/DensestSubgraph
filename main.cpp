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
    Args args;
    args.argsParse(argc, argv);
    std::string graph_path = args.getOption("-g");
    std::string graph_type = args.getOption("-d");
    std::string accuracy = args.getOption("-a");
    std::string epsilon = args.getOption("-e");
    std::string rec_type = args.getOption("-rec");
    std::string alloc_type = args.getOption("-alloc");
    std::string ext_type = args.getOption("-ext");
    std::string ver_type = args.getOption("-ver");

    //合法性检验
    //todo

    Graph graph = Graph(false);
    if (graph_type == "u") {
        graph = Graph(false);
    } else {
        graph = Graph(true);
    }
    graph.loadGraphFromFile(args.getOption("-g"));

    if (accuracy == "a"){
        //todo
    }
    else{
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
            //todo
        }
    }
        return 0;
};