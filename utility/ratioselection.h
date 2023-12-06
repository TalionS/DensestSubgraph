//
// Created by yy on 12/6/23.
//

#ifndef DENSESTSUBGRAPH_RATIOSELECTION_H
#define DENSESTSUBGRAPH_RATIOSELECTION_H

#include "types.h"
#include <vector>
#include <queue>
#include <cmath>
#include <iostream>

class RatioSelection{
public:
    bool ratioSelection(ui vertices_count,
                        std::pair<double, double> &ratio,
                        bool &is_init,
                        bool is_core_dc,
                        bool is_lp_dc,
                        ui s_size,
                        ui t_size);

private:
    struct NormalCompare {
        bool operator()(const std::pair<double, double>& lhs, const std::pair<double, double>& rhs) const {
            return (long long) (lhs.first * rhs.second) > (long long) (rhs.first * lhs.second);
        }
    };


    std::priority_queue<std::pair<double, double>, std::vector<std::pair<double, double>>, NormalCompare> normal_ratio_set_;
    std::priority_queue<double, std::vector<double>, std::greater<>> core_dc_ratio_set_;
private:
    void ratioSetInitialization(ui vertices_count, bool is_core_dc, bool is_lp_dc);
};

#endif //DENSESTSUBGRAPH_RATIOSELECTION_H
