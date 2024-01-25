#include "lp.h"

LinearProgramming::LinearProgramming(bool is_directed, ui type, ui n, ui m, ui sort)
        : // 0 stands for FW 1 stands for FISTA
        is_directed_(is_directed),
        nodes_count_(n),
        edges_count_(m),
        type_(type),
        sort_type(sort) {
    if (is_directed_) {
        r = new std::vector<double>[2];
        alpha.resize(static_cast<unsigned long>(m)), 0;
        for (int i = 0; i < 2; i++) {
            r[i].resize(static_cast<unsigned long>(n), 0);
        }
    } else {
        r = new std::vector<double>[1];
        r[0].resize(static_cast<unsigned long>(n));
        alpha.resize(static_cast<unsigned long>(m));

        if (type_ == 1) {
            beta.resize(static_cast<unsigned long>(m));
        }
    }
}

LinearProgramming::~LinearProgramming() {
    delete[] r;
}

void LinearProgramming::Init(Graph &graph, double ratio) {
    ui n = graph.getVerticesCount();
    ui m = graph.getEdgesCount();
    nodes_count_ = n;
    edges_count_ = m;

    if (is_directed_) {
        ui cnt = 0;
        for (ui i = 0; i < 2; i++)
            r[i].resize(n, 0);
        alpha.resize(m);
        for (VertexID u = 0; u < n; u++) {
            for (auto &v: graph.getOutNeighbors(u)) {
                alpha[cnt].weight_first += 0.5;
                alpha[cnt].weight_second += 0.5;
                alpha[cnt].id_first = u;
                alpha[cnt].id_second = v;
                cnt++;
            }
        }
        for (ui u = 0; u < m; u++) {
            r[0][alpha[u].id_first] += 2 * sqrt(ratio) * alpha[u].weight_first;
            r[1][alpha[u].id_second] += 2 / sqrt(ratio) * alpha[u].weight_second;
        }
        sort(graph);
    } else {
        std::cout << "haha" << std::endl;
        ui cnt = 0;
        r[0].resize(n);
        alpha.resize(m);
        std::cout << "haha" << std::endl;
        for (ui u = 0; u < n; u++) {
            r[0][u] = 0;
            for (auto &v: graph.getNeighbors(u)) {
                if (v > u) continue;
                alpha[cnt].weight_first = 0.5;
                alpha[cnt].weight_second = 0.5;
                alpha[cnt].id_first = u;
                alpha[cnt].id_second = v;
                cnt++;
            }
        }
        std::cout << "haha" << std::endl;
        for (ui u = 0; u < m; u++) {
            r[0][alpha[u].id_first] += alpha[u].weight_first;

            r[0][alpha[u].id_second] += alpha[u].weight_second;
        }
        std::cout << "haha" << std::endl;
        if (type_ == 1) {
            beta.resize(m);
            for (ui i = 0; i < m; i++) {
                beta[i] = alpha[i];
            }
        }
    }
}

void LinearProgramming::Iterate(double learning_rate, double ratio, bool is_synchronous) {
    if (is_directed_) {
        if (is_synchronous) {
            for (ui i = 0; i < nodes_count_; i++) {
                r[0][i] *= (1 - learning_rate);
                r[1][i] *= (1 - learning_rate);
            }
            for (ui i = 0; i < edges_count_; i++) {
                alpha[i].weight_first *= (1 - learning_rate);
                alpha[i].weight_second *= (1 - learning_rate);
                if (r[0][alpha[i].id_first] < r[1][alpha[i].id_second]) {
                    alpha[i].weight_first += learning_rate;
                    r[0][alpha[i].id_first] += 2 * sqrt(ratio) * learning_rate;
                } else {
                    alpha[i].weight_second += learning_rate;
                    r[1][alpha[i].id_second] += 2 / sqrt(ratio) * learning_rate;
                }
            }
        } else {
            std::vector<Alpha> alpha_hat(edges_count_);
            for (ui i = 0; i < edges_count_; i++) {
                alpha[i].weight_first *= (1 - learning_rate);
                alpha[i].weight_second *= (1 - learning_rate);
                if (r[0][alpha[i].id_first] < r[1][alpha[i].id_second]) {
                    alpha[i].weight_first += learning_rate;
//                    r[0][alpha[i].id_first] += 2 * sqrt(ratio) * learning_rate;
                } else {
                    alpha[i].weight_second += learning_rate;
//                    r[1][alpha[i].id_second] += 2 / sqrt(ratio) * learning_rate;
                }
            }
            r[0].assign(nodes_count_, 0);
            r[1].assign(nodes_count_, 0);

//            for (ui u = 0; u < nodes_count_; u++) {
//                r[0][u] = 0;
//                r[1][u] = 0;
//            }
//            for (ui u = 0; u < edges_count_; u++) {
//                r[0][alpha[u].id_first] = 0;
//                r[1][alpha[u].id_second] = 0;
//            }
            for (ui u = 0; u < edges_count_; u++) {
                r[0][alpha[u].id_first] += 2 * sqrt(ratio) * alpha[u].weight_first;
                r[1][alpha[u].id_second] += 2 / sqrt(ratio) * alpha[u].weight_second;
            }
        }
    } else {
        std::vector<Alpha> alpha_hat;
        alpha_hat.resize(edges_count_);
        for (ui i = 0; i < nodes_count_; i++) {
            r[0][i] = (1 - learning_rate) * r[0][i];
        }
//
        for (ui i = 0; i < edges_count_; i++) {
            if (r[0][alpha[i].id_first] < r[0][alpha[i].id_second])
                alpha_hat[i].weight_first = 1;
            else
                alpha_hat[i].weight_second = 1;
        }
        for (ui i = 0; i < edges_count_; i++) {

            alpha[i].weight_first = (1 - learning_rate) * alpha[i].weight_first
                                    + learning_rate * alpha_hat[i].weight_first;
            alpha[i].weight_second = (1 - learning_rate) * alpha[i].weight_second
                                     + learning_rate * alpha_hat[i].weight_second;
        }
        for (ui i = 0; i < nodes_count_; i++) {
            r[0][i] = 0;
        }
        for (ui i = 0; i < edges_count_; i++) {
            r[0][alpha[i].id_first] += alpha[i].weight_first;
            r[0][alpha[i].id_second] += alpha[i].weight_second;
        }
    }
}

void LinearProgramming::FistaIterate(double learning_rate, double t, double ratio) {
    if (is_directed_) {
        //todo
    } else {
        std::vector<Alpha> alpha_new;
        for (ui i = 0; i < edges_count_; i++) {
            beta[i].weight_first = beta[i].weight_first - 2 * learning_rate * r[0][beta[i].id_first];
            beta[i].weight_second = beta[i].weight_second - 2 * learning_rate * r[0][beta[i].id_second];
            if (abs(beta[i].weight_first - beta[i].weight_second) < 1) {
                beta[i].weight_first = (beta[i].weight_first - beta[i].weight_second + 1) / 2;
                beta[i].weight_second = 1 - beta[i].weight_first;
            } else if (beta[i].weight_first - beta[i].weight_second > 0) {
                beta[i].weight_first = 1;
                beta[i].weight_second = 0;
            } else {
                beta[i].weight_first = 0;
                beta[i].weight_second = 1;
            }
        }
        alpha_new = beta;
        for (ui i = 0; i < edges_count_; i++) {
            beta[i].weight_first =
                    alpha_new[i].weight_first + (alpha_new[i].weight_first - alpha[i].weight_first) * (t - 1) / (t + 2);
            beta[i].weight_second = alpha_new[i].weight_second +
                                    (alpha_new[i].weight_second - alpha[i].weight_second) * (t - 1) / (t + 2);
        }
        alpha = alpha_new;
        for (ui i = 0; i < nodes_count_; i++) {
            r[0][i] = 0;
        }
        for (ui i = 0; i < edges_count_; i++) {
            r[0][alpha[i].id_first] += alpha[i].weight_first;
            r[0][alpha[i].id_second] += alpha[i].weight_second;
        }
    }
}

void LinearProgramming::sort(Graph &graph) {
    if (is_directed_) {
        if (sort_type == 1) {
            std::vector<std::vector<ui>> degrees(2);
            degrees[0] = graph.getOutDegrees();
            degrees[1] = graph.getInDegrees();
            std::sort(alpha.begin(), alpha.end(), [this, &degrees](Alpha a, Alpha b) -> bool {
                return degrees[0][a.id_first] * degrees[1][a.id_second] <
                       degrees[0][b.id_first] * degrees[1][b.id_second];
            });
        } else if (sort_type == 2) {
            std::vector<std::vector<ui>> degrees(2);
            degrees[0] = graph.getOutDegrees();
            degrees[1] = graph.getInDegrees();
            std::sort(alpha.begin(), alpha.end(), [this, &degrees](Alpha a, Alpha b) -> bool {
                return degrees[0][a.id_first] + degrees[1][a.id_second] <
                       degrees[0][b.id_first] + degrees[1][b.id_second];
            });
        }
    } else {

    }
}