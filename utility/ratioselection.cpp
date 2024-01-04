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
}

bool RatioSelection::ratioSelection(ui vertices_count,
                                    std::pair<double, double> &ratio,
                                    bool &is_init,
                                    bool is_core_dc,
                                    bool is_lp_dc,
                                    ui s_size,
                                    ui t_size,
                                    double ratio_o,
                                    double ratio_p){
    if (!is_init) {
        is_init = true;
        ratioSetInitialization(vertices_count, is_core_dc, is_lp_dc);
        if(is_core_dc || is_lp_dc){
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
        while(!dc_ratio_set_.empty() || !is_pushed){
            if(!is_pushed) {
                is_pushed = true;
                double mid = (ratio.first + ratio.second) / 2;
                double b = mid;
                if (s_size && t_size)
                    b = (double) s_size / t_size;
//                else
//                    continue;
//                double c = mid * mid / b;
//                if (b > c)
//                    std::swap(b, c);
//                if (ratio.first < b && b < ratio.second){
//                    dc_ratio_set_.push(ratio.first);
//                    dc_ratio_set_.push(b);
//                }
//                if(ratio.second > c && c > ratio.first){
//                    dc_ratio_set_.push(ratio.second);
//                    dc_ratio_set_.push(c);
//                }
//                printf("%d, %d\n", s_size, t_size);
                double c = sqrt(b) / sqrt(mid) + sqrt(mid) / sqrt(b);
                if(b > mid){
                    double mid_cover = (-2 * mid + mid * c * c - mid * sqrt(pow(c, 4) - 4 * c * c)) / 2;
                    if (b < mid_cover)
                        std::swap(b ,mid_cover);
                    if (ratio.first < mid_cover && mid_cover < ratio.second) {
                        dc_ratio_set_.push(ratio.first);
                        dc_ratio_set_.push(mid_cover);
                    }
                    if (ratio.first <  b && b < ratio.second){
                        dc_ratio_set_.push(b);
                        dc_ratio_set_.push(ratio.second);
                    }
                }
                else{
                    double mid_cover = (-2 * mid + mid * c * c + mid * sqrt(pow(c, 4) - 4 * c * c)) / 2;
                    if (b > mid_cover)
                        std::swap(b, mid_cover);
                    if (ratio.first < mid_cover && mid_cover < ratio.second){
                        dc_ratio_set_.push(ratio.second);
                        dc_ratio_set_.push(mid_cover);
                    }
                    if (ratio.second > b && b > ratio.first){
                        dc_ratio_set_.push(b);
                        dc_ratio_set_.push(ratio.first);
                    }
                }
//                printf("(%f, %f, %f, %f)", ratio.first, b, c, ratio.second);
//                if (ratio.first <= b){
//                    dc_ratio_set_.push(b);
//                    dc_ratio_set_.push(ratio.first);
//                }
//                if (ratio.second >= c){
//                    dc_ratio_set_.push(c);
//                    dc_ratio_set_.push(ratio.second);
//                }
//                dc_ratio_set_.push(b);
//                dc_ratio_set_.push(c);
            }
            if(dc_ratio_set_.empty())
                continue;
            ratio.first = dc_ratio_set_.top();
            dc_ratio_set_.pop();
            ratio.second = dc_ratio_set_.top();
            dc_ratio_set_.pop();
            if (ratio.first + 1.0 / vertices_count > ratio.second)
                continue;
            return true;
        }
    }
    if(is_lp_dc){
        bool is_pushed = false;
        while (!dc_ratio_set_.empty() || !is_pushed){
            if(!is_pushed){
                is_pushed = true;
                if (ratio.second > ratio_o && ratio_o > ratio.first){
                    dc_ratio_set_.push(ratio_o);
                    dc_ratio_set_.push(ratio.first);
                }
                if(ratio.first < ratio_p && ratio_p < ratio.second){
                    dc_ratio_set_.push(ratio_p);
                    dc_ratio_set_.push(ratio.second);
                }
            }
            ratio.first = dc_ratio_set_.top();
            dc_ratio_set_.pop();
            ratio.second = dc_ratio_set_.top();
            dc_ratio_set_.pop();
            if (ratio.first + 1.0 / vertices_count > ratio.second)
                continue;
            return true;
        }
    }
    return false;
}
