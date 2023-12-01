//
// Created by yy on 11/30/23.
//

#ifndef DENSESTSUBGRAPH_ARGS_H
#define DENSESTSUBGRAPH_ARGS_H

#include <map>
#include <vector>

class Args {
private:
    std::map<std::string, std::string> args_;
    std::vector<std::string> options_ = {
            "-g",       //graph file path
            "-d",       //type of graph, u for undirected, d for directed
            "-a",       //accuracy, e for exact, a for approximate
            "-e",       //epsilon
            "-rec",     //type of reduction
            "-alloc",   //type of allocation
            "-ext",     //type of extraction
            "-ver"      //type of verification

    };

public:
    Args();

    void argsParse(int argc, char **argv);

    std::string getOption(const std::string &option);
};



//        flow_exact,
//        core_exact,
//        lp_exact


//        greedy,
//        core_app,
//        greedtpp,
//        fista,
//        flow_app


//        flow_exact,
//        dc_exact,
//        cp_exact
//

//        bs_app,
//        ks_app,
//        core_app,
//        w_core_app,
//        vw_app,
//        cp_app

#endif //DENSESTSUBGRAPH_ARGS_H
