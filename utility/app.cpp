#include "app.h"

struct rank{
    ui x,y;
};
bool cmp(rank a, rank b){
    return a.y > b.y;
};
CoreApp::CoreApp(ui nodes_count, ui pos) : // 0 stands for FW 1 stands for FISTA
        nodes_count(nodes_count),
        pos(pos) {
    k = 0;
    opt = 0;
    w.resize(nodes_count);
}

void CoreApp::Init(Graph &graph) {
    deg = graph.CoreOrder(graph);
    nodes_count = graph.getVerticesCount();
    id.resize(nodes_count);
    w.resize(nodes_count);
    selected.resize(nodes_count);
    new_edge = new std::vector<ui>[nodes_count];
    for(ui i = 0; i < nodes_count; i++){
        new_edge[i].clear();
    }
    size = 100;
    /*std::vector<ui> deg = graph.getDegrees();
    std::vector<ui> id,rev;
    id.resize(nodes_count);
    rev.resize(nodes_count);
    w.resize(nodes_count);
    std::vector<rank> tmp;
    tmp.resize(nodes_count);
    for(ui i = 0; i < nodes_count; i++){
        tmp[i].x = i;
        tmp[i].y = deg[i];
    }
    //std::sort(tmp.begin(), tmp.end(), cmp);
    for(ui i = 0; i < nodes_count; i++){
        id[i] = tmp[i].x;
        rev[tmp[i].x] = i;
    }
    while(pos + 1 < nodes_count && deg[id[pos]] == deg[id[pos + 1]]) pos++;
    */
}

PKMC::PKMC(){ // 0 stands for FW 1 stands for FISTA
    hmax = 0;
    s = 0;
}
void PKMC::Init(Graph &graph){
    ui vertices_count = graph.getVerticesCount();
    h = graph.getDegrees();
    hmax = 0;
    for(ui i = 0; i < vertices_count; i++){
        hmax = std::max(hmax, h[i]);
    }
    s = 0;
    for(ui i = 0; i < vertices_count; i++){
        s += (h[i] == hmax);
    }
}

void Greedy::Init(Graph &graph){
    ui vertices_count = graph.getVerticesCount();
    h.resize(vertices_count);
    for(ui i = 0; i < vertices_count; i++){
        h[i] = 0;
    }
}