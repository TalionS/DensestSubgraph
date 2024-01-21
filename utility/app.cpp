#include "app.h"

struct rank{
    ui x,y;
}
bool cmp(rank a, rank b){
    return a.y > b.y;
}
CoreApp::CoreApp(ui nodes_count, ui pos) : // 0 stands for FW 1 stands for FISTA
        nodes_count(nodes_count),
        pos(pos) {
    id.resize(nodes_count);
    w.resize(nodes_count);
}

void CoreApp::Init(Graph &graph) {
    nodes_count = graph.getVerticesCount();
    vector<ui> deg = graph.getDegrees();
    id.resize(nodes_count);
    w.resize(nodes_count);
    vector<rank> tmp;
    tmp.resize(nodes_count);
    for(ui i = 0; i < nodes_count; i++){
        tmp[i].x = i;
        tmp[i].y = deg[i];
    }
    sort(tmp.begin(), tmp.end(), cmp);
    for(ui i = 0; i < nodes_count; i++){
        id[i] = tmp[i].x;
    }
    while(pos + 1 < nodes_count && deg[id[pos]] == deg[id[pos + 1]]) pos++;
}