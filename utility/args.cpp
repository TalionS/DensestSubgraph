//
// Created by yy on 11/30/23.
//

#include "args.h"
#include <vector>
#include <string>
#include <algorithm>

Args::Args() {

    for (const auto &opt: options_) {
        args_[opt] = "default";
    }
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
        args_[argv[i++]] = argv[i];
    }
}

std::string Args::getOption(const std::string &option) {
    if (std::find(options_.begin(), options_.end(), option) == options_.end()) {
        //error
        exit(-1);
    }
    return args_[option];
}
