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
                                    bool is_dc, double ratio_o, double ratio_p, double density, double epsilon) {
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
            ratio.first = std::max(ratio.first, (density / 2 / max_degree[0]) * (density / 2 / max_degree[0]));
            ratio.second = std::min(ratio.second, (2 * max_degree[1] / density) * (2 * max_degree[1] / density));
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
//    if (is_core_dc) {
//        bool is_pushed = false;
//        while (!dc_ratio_set_.empty() || !is_pushed) {
//            if (!is_pushed) {
//                is_pushed = true;
//                double mid = (ratio.first + ratio.second) / 2;
//                double b = mid;
//                if (s_size && t_size)
//                    b = (double) s_size / t_size;
//                double c = sqrt(b) / sqrt(mid) + sqrt(mid) / sqrt(b);
//                if (b > mid) {
//                    double mid_cover = (-2 * mid + mid * c * c - mid * sqrt(pow(c, 4) - 4 * c * c)) / 2;
//                    if (b < mid_cover)
//                        std::swap(b, mid_cover);
//                    if (ratio.first < mid_cover && mid_cover < ratio.second) {
//                        dc_ratio_set_.push(ratio.first);
//                        dc_ratio_set_.push(mid_cover);
//                    }
//                    if (ratio.first < b && b < ratio.second) {
//                        dc_ratio_set_.push(b);
//                        dc_ratio_set_.push(ratio.second);
//                    }
//                } else {
//                    double mid_cover = (-2 * mid + mid * c * c + mid * sqrt(pow(c, 4) - 4 * c * c)) / 2;
//                    if (b > mid_cover)
//                        std::swap(b, mid_cover);
//                    if (ratio.first < mid_cover && mid_cover < ratio.second) {
//                        dc_ratio_set_.push(ratio.second);
//                        dc_ratio_set_.push(mid_cover);
//                    }
//                    if (ratio.second > b && b > ratio.first) {
//                        dc_ratio_set_.push(b);
//                        dc_ratio_set_.push(ratio.first);
//                    }
//                }
//            }
//            if (dc_ratio_set_.empty())
//                continue;
//            ratio.first = dc_ratio_set_.top();
//            dc_ratio_set_.pop();
//            ratio.second = dc_ratio_set_.top();
//            dc_ratio_set_.pop();
//            ratio.first = std::max(ratio.first, (density / 2 / max_degree[0]) * (density / 2 / max_degree[0]));
//            ratio.second = std::min(ratio.second, (2 * max_degree[1] / density) * (2 * max_degree[1] / density));
//            if (ratio.first + 1.0 / vertices_count > ratio.second)
//                continue;
//            return true;
//        }
//    }
    if (is_dc) {
        bool is_pushed = false;
        while (!dc_ratio_set_.empty() || !is_pushed) {
            if (!is_pushed) {
                is_pushed = true;
                if (ratio.second > ratio_o && ratio_o > ratio.first) {
                    dc_ratio_set_.push(ratio_o);
                    dc_ratio_set_.push(ratio.first);
                }
                if (ratio.first < ratio_p && ratio_p < ratio.second) {
                    dc_ratio_set_.push(ratio_p);
                    dc_ratio_set_.push(ratio.second);
                }
            }
            if (dc_ratio_set_.empty())
                continue;
            ratio.first = dc_ratio_set_.top();
            dc_ratio_set_.pop();
            ratio.second = dc_ratio_set_.top();
            dc_ratio_set_.pop();
            ratio.first = std::max(ratio.first, (density / 2 / max_degree[0]) * (density / 2 / max_degree[0]));
            ratio.second = std::min(ratio.second, (2 * max_degree[1] / density) * (2 * max_degree[1] / density));
            if (ratio.first + 1.0 / vertices_count > ratio.second)
                continue;
            return true;
        }
    }
    return false;
}
