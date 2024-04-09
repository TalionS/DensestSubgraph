//
// Created by yy on 12/6/23.
//

#include "ratioselection.h"

RatioSelection::RatioSelection(Graph &graph) {
    std::vector<std::vector<ui>> deg(2);
    deg[0] = graph.getOutDegrees();
    deg[1] = graph.getInDegrees();
    max_degree.resize(2);
    max_degree[0] = *std::max_element(deg[0].begin(), deg[0].end());
    max_degree[1] = *std::max_element(deg[1].begin(), deg[1].end());
}

void RatioSelection::ratioSetInitialization(ui vertices_count, bool is_vw, bool is_dc) {
    if (!is_dc && !is_vw) {
        for (ui i = 1; i <= vertices_count; i++) {
            normal_ratio_set_.push(std::make_pair(1, i));
        }
    }
}

bool RatioSelection::ratioSelection(ui vertices_count, std::pair<double, double> &ratio, bool &is_init, bool is_vw,
                                    bool is_dc, double &ratio_o, double &ratio_p, double density, double epsilon,
                                    bool is_res, ui res_width, bool is_map) {
    if (!is_init) {
        is_init = true;
        ratioSetInitialization(vertices_count, false, is_dc);
        if (is_vw) {
            ratio.first = 1;
            ratio.second = vertices_count;
            return true;
        }
        if (is_dc) {
            ratio.first = 1.0 / vertices_count;
            ratio.second = vertices_count;
            ratio.first = std::max(ratio.first, (density / max_degree[0]) * (density / max_degree[0]));
            ratio.second = std::min(ratio.second, (max_degree[1] / density) * (max_degree[1] / density));
            return true;
        }
    }
    if (!is_dc && !is_vw) {
        while (!normal_ratio_set_.empty()) {
            ratio = normal_ratio_set_.top();
            normal_ratio_set_.pop();
            if (ratio.first > vertices_count)
                continue;
            normal_ratio_set_.push(std::make_pair(ratio.first + 1, ratio.second));
            if ((long long) (ratio.first * normal_ratio_set_.top().second) ==
                (long long) (ratio.second * normal_ratio_set_.top().first))
                continue;
            return true;
        }
    }
    if (is_vw) {
        ratio.first *= (2 + 2 * epsilon) / (2 + epsilon);
        ratio.first *= (2 + 2 * epsilon) / (2 + epsilon);
        if (ratio.first / ratio.second > vertices_count)
            return false;
        return true;
    }
    if (is_dc) {
        bool is_pushed = false;
        double c, res_ratio_l, res_ratio_r;
        if (is_map) {
            if (ratio.first < 1 && ratio.second > 1) {
                c = 1;
            } else if (ratio.second <= 1) {
                c = (ratio.first + ratio.second) / 2;
            } else if (ratio.first >= 1) {
                c = 2 / (1 / ratio.first + 1 / ratio.second);
            }
        } else
            c = (ratio.first + ratio.second) / 2;
        if (is_res) {
            res_ratio_l = std::max(ratio.first, c / res_width);
            res_ratio_r = std::min(ratio.second, c * res_width);
        }
        while (!dc_ratio_set_.empty() || !is_pushed) {
            if (!is_pushed) {
                is_pushed = true;
                if (ratio_o != 0 && is_res) {
                    ratio_o = std::max(ratio_o, res_ratio_l);
                    ratio_p = std::min(ratio_p, res_ratio_r);
                }
                if (ratio_o == 0) {
//                    ratio_o = (ratio.first + ratio.second) / 2;
//                    ratio_p = ratio_o;
                    ratio_o = c;
                    ratio_p = c;
                }

                if (ratio.second > ratio_o && ratio_o > ratio.first) {
                    dc_ratio_set_.push(ratio_o);
                    dc_ratio_set_.push(ratio.first);
                }
                if (ratio.first < ratio_p && ratio_p < ratio.second) {
                    dc_ratio_set_.push(ratio_p);
                    dc_ratio_set_.push(ratio.second);
                }
                ratio_o = 0;
                ratio_p = 0;
            }
            if (dc_ratio_set_.empty())
                continue;
            ratio.first = dc_ratio_set_.top();
            dc_ratio_set_.pop();
            ratio.second = dc_ratio_set_.top();
            dc_ratio_set_.pop();
            ratio.first = std::max(ratio.first, (density / max_degree[0]) * (density / max_degree[0]));
            ratio.second = std::min(ratio.second, (max_degree[1] / density) * (max_degree[1] / density));
            if (ratio.first + 1.0 / vertices_count > ratio.second)
                continue;
            return true;
        }
    }
    return false;
}
