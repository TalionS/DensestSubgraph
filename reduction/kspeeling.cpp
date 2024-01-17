#include "kspeeling.h"

KSPeeling::KSPeeling(Graph &graph) {
    vert_.resize(2);
    bin_.resize(2);
    pos_.resize(2);
    degrees_.resize(2);
    auto vertices_count = graph.getVerticesCount();
    degrees_[0] = graph.getOutDegrees();
    degrees_[1] = graph.getInDegrees();
    for (int i = 0; i < 2; i++) {
        vert_[i].clear();
        bin_[i].clear();
        pos_[i].clear();

        vert_[i].resize(vertices_count, 0);
        bin_[i].resize(vertices_count + 1, 0);
        pos_[i].resize(vertices_count, 0);
        for (VertexID v = 0; v < vertices_count; v++) {
            ++bin_[i][degrees_[i][v]];
        }
        ui start = 0;
        for (int d = 0; d <= vertices_count; d++) {
            ui num = bin_[i][d];
            bin_[i][d] = start;
            start += num;
        }
        for (VertexID v = 0; v < vertices_count; v++) {
            pos_[i][v] = bin_[i][degrees_[i][v]]++;
            vert_[i][pos_[i][v]] = v;
        }
        for (int d = vertices_count; d > 0; d--) {
            bin_[i][d] = bin_[i][d - 1];
        }
        bin_[i][0] = 0;
    }
}
using Heap = boost::heap::fibonacci_heap<std::pair<ui, VertexID>>;
void KSPeeling::peeling(double &density, std::vector<std::vector<VertexID>> &vertices) {

}

inline void KSPeeling::decDeg(ui cur, VertexID t, std::vector<std::vector<VertexID>> &vert_copy, std::vector<std::vector<ui>> &bin_copy, std::vector<std::vector<ui>> &pos_copy,
                           std::vector<std::vector<ui>> &degrees_copy) {
    ui dt = degrees_copy[cur][t], pt = pos_copy[cur][t];
    ui pw = bin_copy[cur][dt];
    VertexID w = vert_copy[cur][pw];
    if (t != w) {
        pos_copy[cur][t] = pw; vert_copy[cur][pt] = w;
        pos_copy[cur][w] = pt; vert_copy[cur][pw] = t;
    }
    ++bin_copy[cur][dt];
    --degrees_copy[cur][t];
}
