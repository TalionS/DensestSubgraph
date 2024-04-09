#ifndef DENSESTSUBGRAPH_HEAP_H
#define DENSESTSUBGRAPH_HEAP_H

#include <vector>
#include <iostream>
#include "types.h"
#include "graph.h"
#include <fstream>
#include <string>
#include <algorithm>
#include <random>
#include <cmath>
#include <cassert>


struct HeapElement{
    ui id;
    double val;
};

class BinaryHeap {
private:
    std::vector<HeapElement> heap;
    std::vector<ui> pos;
    void heapifyUp(ui index) {
        while (index > 1) {
            ui parent = index / 2;
            if (heap[parent].val < heap[index].val)
                break;
            std::swap(heap[parent], heap[index]);
            pos[heap[parent].id] = parent;
            pos[heap[index].id] = index;
            index = parent;
        }
    }

    void heapifyDown(ui index) {
        ui smallest = index;
        ui left = 2 * index;
        ui right = 2 * index + 1;
        if (left < heap.size() && heap[left].val < heap[smallest].val)
            smallest = left;

        if (right < heap.size() && heap[right].val < heap[smallest].val)
            smallest = right;

        if (smallest != index) {
            std::swap(heap[index], heap[smallest]);
            pos[heap[index].id] = index;
            pos[heap[smallest].id] = smallest;
            heapifyDown(smallest);
        }
    }

public:
    BinaryHeap() {}

    void resize(ui nodes_count){
        heap.resize(nodes_count + 1);
        pos.resize(nodes_count);
    }

    void init(std::vector<double> &h){
        for(ui i = 0; i < h.size(); i++){
            heap[i + 1].id = i;
            heap[i + 1].val = h[i];
            pos[i] = i + 1;
            heapifyUp(i + 1);
        }
    }

    HeapElement pop() {
        HeapElement res = heap[1];
        ui sz = heap.size() - 1;
        std::swap(heap[1], heap[sz]);
        pos[heap[1].id] = 1;
        pos[heap[sz].id] = sz;
        heap.pop_back();
        heapifyDown(1);
        return res;
    }

    void modify(ui id){
        ui index = pos[id];
        heap[index].val--;
        heapifyUp(index);
    }

    void dec(ui id, double delt, double minval){
        ui index = pos[id];
        heap[index].val = std::max(heap[index].val - delt, minval);
        heapifyUp(index);
    }
};
#endif //DENSESTSUBGRAPH_HEAP_H