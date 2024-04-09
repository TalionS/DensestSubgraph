//
// Created by yy on 12/1/23.
//

#ifndef DENSESTSUBGRAPH_ALLOCATION_H
#define DENSESTSUBGRAPH_ALLOCATION_H

#include "../utility/flownetwork.h"
#include "../utility/graph.h"
#include "../utility/xycore.h"
#include "../utility/lp.h"
#include "../utility/app.h"
#include "../utility/heap.h"
#include "wcore.h"
#include <algorithm>
#include <cmath>
#include <boost/heap/fibonacci_heap.hpp>
using Heap = boost::heap::fibonacci_heap<std::pair<int, VertexID>>;

//[this] (std::pair<ui, VertexID> a, std::pair<ui, VertexID> b)->bool {
//return a.first < b.first || (a.first == b.first && a.first <= b.first);}
class Allocation{
public:
    void flowExactAllocation(Graph &graph, FlowNetwork &flow, std::pair<double, double> ratio, double l, double r,
                             bool is_dc, bool is_map);
    void xyCoreApproAllocation(Graph &graph, std::pair<ui, ui> &max_core_num_pair);
    void wCoreApproAllocation(Graph &graph, WCore &w_core, std::pair<ui, ui> &max_core_num_pair);
    void directedKSApproAllocation(Graph &graph, std::vector<Heap> &heap,
                                   std::vector<std::vector<Heap::handle_type>> &handles,
                                   std::vector<std::vector<bool>> &is_peeled,
                                   ui &edges_count,
                                   bool &is_init);
    void directedFixedKSApproAllocation(Graph &graph, std::pair<double, double> ratio, ui cur, std::vector<Heap> &heap,
                                        std::vector<std::vector<Heap::handle_type>> &handles,
                                        std::vector<std::vector<bool>> &is_peeled,
                                        ui &edges_count);
    void directedBSApproAllocation(Graph &graph, std::pair<double, double> ratio, std::vector<Heap> &heap,
                                   std::vector<std::vector<Heap::handle_type>> &handles,
                                   std::vector<std::vector<bool>> &is_peeled,
                                   ui &edges_count,
                                   bool &is_init);
    void directedPMApproAllocation(Graph &graph, std::pair<double, double> ratio, double epsilon, ui &edges_count,
                                   std::vector<std::vector<VertexID>> &vertices, std::vector<std::vector<ui>> &degrees,
                                   bool &is_init);
    void
    directedCPAllocation(Graph &graph, LinearProgramming &lp, ui &iter_num, bool &is_init,
                         std::pair<double, double> ratios, bool is_synchronous, bool is_exp, bool is_map);
//    void directedVWApproAllocation(Graph &graph, LinearProgramming &lp, ui T, bool &is_init, std::pair<double, double> ratios)
    void UndirectedflowExactAllocation(Graph &graph, FlowNetwork &flow, double l, double r);
    void UndirectedlpAllocation(Graph &graph, LinearProgramming &lp, ui T, bool is_synchronous = false);
    void UndirectedFistaAllocation(Graph &graph, LinearProgramming &lp, ui T, bool is_synchronous = false);
    void UndirectedMWUAllocation(Graph &graph, LinearProgramming &lp, ui T, bool is_synchronous = false);
    void UndirectedCoreAppAllocation(Graph &graph, CoreApp &ca);
    void UndirectedGreedyAllocation(Graph &graph);
    void UndirectedGreedyppAllocation(Graph &graph, Greedy &gr);
    void UndirectedFlowAppAllocation(Graph &graph);
private:
    XYCore xycore;
};

#endif //DENSESTSUBGRAPH_ALLOCATION_H
