//
// Created by yy on 12/1/23.
//

#include "verification.h"

bool Verification::flowExactVerification(Graph &graph, double l, double r) {
    ui n = graph.getVerticesCount();
    double bias = 1.0 / sqrt((double) n * (n - 1)) - 1.0 / n;

    if (r - l > bias) return true;
    else return false;
}