#include "lp.h"

LinearProgramming::LinearProgramming(bool is_directed, ui n, ui m) :
        is_directed_(is_directed),
        nodes_count_(n),
        edges_count_(m) {
    if (is_directed_) {
        r = new std::vector<double>[2];
        alpha = new std::vector<Alpha>[1];
        alpha[0].resize(static_cast<unsigned long>(m)), 0;
        for (int i = 0; i < 2; i++) {
            r[i].resize(static_cast<unsigned long>(n), 0);
        }
    } else {
        r = new std::vector<double>[1];
        alpha = new std::vector<Alpha>[1];
        r[0].resize(static_cast<unsigned long>(n));
        alpha[0].resize(static_cast<unsigned long>(m));
    }
}

LinearProgramming::~LinearProgramming() {
    delete[] r;
    delete[] alpha;
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
        alpha[0].resize(m);
        for (VertexID u = 0; u < n; u++) {
            for (auto &v: graph.getOutNeighbors(u)) {
                alpha[0][cnt].weight_first += 0.5;
                alpha[0][cnt].weight_second += 0.5;
                alpha[0][cnt].id_first = u;
                alpha[0][cnt].id_second = v;
                cnt++;
            }
        }
        for (ui u = 0; u < m; u++) {
            r[0][alpha[0][u].id_first] += 2 * sqrt(ratio) * alpha[0][u].weight_first;
            r[1][alpha[0][u].id_second] += 2 / sqrt(ratio) * alpha[0][u].weight_second;
        }
//        printf("alpha\n");
//        for(auto a: alpha[0])
//            std::cout << "(" << a.id_first << a.id_second << ")" << a.weight_first << " " << a.weight_second <<std::endl;
//        printf("r\n");
//        for (auto r_: r[0])
//            std::cout << r_ << " ";
//        std::cout << std::endl;
//        for (auto r_ : r[1])
//            std::cout << r_ << " ";
//        std::cout << std::endl;
//        printf("Init finished.\n");
    } else {
        ui cnt = 0;
        r[0].resize(n);
        alpha[0].resize(m);
        for (ui u = 0; u < n; u++) {
            r[0][u] = 0;
            for (auto &v: graph.getNeighbors(u)) {
                if (v > u) continue;
                alpha[0][cnt].weight_first = 0.5;
                alpha[0][cnt].weight_second = 0.5;
                alpha[0][cnt].id_first = u;
                alpha[0][cnt].id_second = v;
                cnt++;
            }
        }
        for (ui u = 0; u < m; u++) {
            r[0][alpha[0][u].id_first] += alpha[0][u].weight_first;
            r[1][alpha[0][u].id_second] += alpha[0][u].weight_second;
        }
    }
}

void LinearProgramming::Iterate(double learning_rate, double ratio) {
    if (is_directed_) {
        for (ui i = 0; i < nodes_count_; i++) {
            r[0][i] *= (1 - learning_rate);
            r[1][i] *= (1 - learning_rate);
        }
        random_shuffle(alpha[0].begin(), alpha[0].end());
        for (ui i = 0; i < edges_count_; i++) {
            alpha[0][i].weight_first *= (1 - learning_rate);
            alpha[0][i].weight_second *= (1 - learning_rate);
            if (r[0][alpha[0][i].id_first] < r[1][alpha[0][i].id_second]) {
                alpha[0][i].weight_first += learning_rate;
                r[0][alpha[0][i].id_first] += 2 * sqrt(ratio) * learning_rate;
            } else if (r[0][alpha[0][i].id_first] > r[1][alpha[0][i].id_second]) {
                alpha[0][i].weight_second += learning_rate;
                r[1][alpha[0][i].id_second] += 2 / sqrt(ratio) * learning_rate;
            } else {
                if (ratio < 1) {
                    alpha[0][i].weight_first += learning_rate;
                    r[0][alpha[0][i].id_first] += 2 * sqrt(ratio) * learning_rate;
                } else {
                    alpha[0][i].weight_second += learning_rate;
                    r[1][alpha[0][i].id_second] += 2 / sqrt(ratio) * learning_rate;
                }
            }
        }
    } else {
        std::vector <Alpha> alpha_hat;
        alpha_hat.resize(edges_count_);
        for (ui i = 0; i < nodes_count_; i++) {
            r[0][i] = (1 - learning_rate) * r[0][i];
        }
        random_shuffle(alpha[0].begin(), alpha[0].end());
        for (ui i = 0; i < edges_count_; i++) {
            alpha[0][i].weight_first = (1 - learning_rate) * alpha[0][i].weight_first;
            alpha[0][i].weight_second = (1 - learning_rate) * alpha[0][i].weight_second;
            if (r[0][alpha[0][i].id_first] < r[0][alpha[0][i].id_second]) {
                alpha[0][i].weight_first += learning_rate;
                r[0][alpha[0][i].id_first] += learning_rate;
            } else {
                alpha[0][i].weight_second += learning_rate;
                r[0][alpha[0][i].id_second] += learning_rate;
            }
        }
        /*
        for(ui i = 0; i < edges_count_; i++){
            if(r[0][alpha[0][i].id_first] < r[0][alpha[0][i].id_second])
                alpha_hat[i].weight_first = 1;
            else
                alpha_hat[i].weight_second = 1;
        }
        for(ui i = 0; i < edges_count_; i++){
            alpha[0][i].weight_first = (1 - learning_rate) * alpha[0][i].weight_first 
                + learning_rate * alpha_hat[i].weight_first;
            alpha[0][i].weight_second = (1 - learning_rate) * alpha[0][i].weight_second 
                + learning_rate * alpha_hat[i].weight_second;
        }
        for(ui i = 0; i < nodes_count_; i++){
            r[0][i] = 0;
        }
        for(ui i = 0; i < edges_count_; i++){
            r[0][alpha[0][i].id_first] += alpha[0][i].weight_first;
            r[0][alpha[0][i].id_second] += alpha[0][i].weight_second;
        }
        */
    }
}