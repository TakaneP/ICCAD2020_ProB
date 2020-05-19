#include<bits/stdc++.h>

#ifndef SEGMENTTREE
#define SEGMENTTREE

class RoutingGraph;

struct SegmentTreeNode{
    int lzTag;
    int minValue;
};

class SegmentTree{
public:
    SegmentTree() = delete;
    SegmentTree(RoutingGraph& main): graph(main) {}
    RoutingGraph& graph;
    std::vector<std::vector<std::vector<SegmentTreeNode>>> node; //layer, n = |row| or |column|, 2*2^ceil(logn)
    void build_ini(void);
    void build(int treeNodeIndex, int lowerBound, int upperBound, int layer, int rowOrColIndex);
    void pushup(int treeNodeIndex, int layer, int rowOrColIndex);
    void pushdown(int treeNodeIndex, int layer, int rowOrColIndex);
    int get_remaining_supply(int startIndex, int endIndex, int layer, int rowOrColIndex);
    int query(int treeNodeIndex, int lowerBound, int upperBound, int startIndex, int endIndex, int layer, int rowOrColIndex);
    void update(int treeNodeIndex, int lowerBound, int upperBound, int startIndex, int endIndex, int layer, int rowOrColIndex, int delta);

};



#endif
