//
// Created by yy on 12/6/23.
//

#include "ratioselection.h"

void RatioSelection::ratioSetInitialization(ui vertices_count, bool is_core_dc, bool is_lp_dc) {
    if (!is_core_dc && !is_lp_dc){
        for (ui i = 0; i < vertices_count; i++){
            normal_ratio_set_.push(std::make_pair(1, i));
        }
    }

    if(is_core_dc){
        core_dc_ratio_set_.push(1.0 / vertices_count);
        core_dc_ratio_set_.push(vertices_count);
    }
}

bool RatioSelection::ratioSelection(ui vertices_count,
                                    std::pair<double, double> &ratio,
                                    bool &is_init,
                                    bool is_core_dc,
                                    bool is_lp_dc,
                                    ui s_size,
                                    ui t_size){
    if (!is_init) {
        is_init = true;
        ratioSetInitialization(vertices_count, is_core_dc, is_lp_dc);
        if(is_core_dc){
            ratio.first = 1.0 / vertices_count;
            ratio.second = vertices_count;
            return true;
        }

    }
    if(!is_core_dc && !is_lp_dc){
        while(!normal_ratio_set_.empty()){
            ratio = normal_ratio_set_.top();
            normal_ratio_set_.pop();
            if (ratio.first > vertices_count)
                continue;
            normal_ratio_set_.push(std::make_pair(ratio.first + 1, ratio.second));
            if ((long long) (ratio.first * normal_ratio_set_.top().second) == (long long) (ratio.second * normal_ratio_set_.top().first))
                continue;
            return true;
        }
    }
    if(is_core_dc){
        bool is_pushed = false;
        while(!core_dc_ratio_set_.empty()){
            if(!is_pushed) {
                is_pushed = true;
                double mid = (ratio.first + ratio.second) / 2;
                double b = mid;
                if (s_size && t_size)
                    b = (double) s_size / t_size;
                double c = sqrt(b) / sqrt(mid) + sqrt(mid) / sqrt(b);
                printf("(%f, %f)", std::min(b, c), std::max(b, c));
                core_dc_ratio_set_.push(b);
                core_dc_ratio_set_.push(c);
            }
            ratio.first = core_dc_ratio_set_.top();
            core_dc_ratio_set_.pop();
            ratio.second = core_dc_ratio_set_.top();
            core_dc_ratio_set_.pop();
            if (ratio.first + 1.0 / vertices_count > ratio.second)
                continue;
            return true;
        }
    }
    return false;
}
