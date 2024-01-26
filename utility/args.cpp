//
// Created by yy on 11/30/23.
//

#include "args.h"


Args::Args() {
    args_[options_[0]] = "./";
    args_["-t"] = "u";
    args_["-a"] = "e";
    args_["-eps"] = "0";
    args_["-red"];
    args_["-alloc"];
    args_["-ext"];
    args_["-ver"];
    args_["-o"] = "0";
    args_["-seq"] = "f";
    args_["-vw"] = "f";
    args_["-lr"] = "0";
    args_["-p"] = "f";
    args_["-exp"] = "t";
    args_["-it"] = "100";
    args_["-dc"] = "t";
}


void Args::argsParse(int argc, char **argv) {
    for (int i = 1; i < argc; i++) {
        if (std::find(options_.begin(), options_.end(), argv[i]) == options_.end()) {
            //error
            exit(-1);
        }
        if (i + 1 > argc) {
            //error
            exit(-1);
        }
        args_[argv[i]] = argv[i + 1];
        i++;
    }
}

std::string Args::getOption(const std::string &option) {
    if (std::find(options_.begin(), options_.end(), option) == options_.end()) {
        //error
        exit(-1);
    }
    return args_[option];
}
