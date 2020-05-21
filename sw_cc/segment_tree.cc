#include<bits/stdc++.h>
#include "data_structure.h"
#include "segment_tree.h"

using namespace std;
//Trigger Build function
void SegmentTree::build_ini(void) {
    for(int i = 0; i < graph.layer; ++i) {
        if(i & 1){ // vertical
            for(int j = 0; j < graph.column; j++) {
                build(1, 0, graph.row-1, i, j);
            }
        }
        else { //horizontal
            for(int j = 0; j < graph.row; j++) {
                build(1, 0, graph.column - 1, i, j);
            }
        }
    }
}
//Recursive build segment tree
void SegmentTree::build(int treeNodeIndex, int lowerBound, int upperBound, int layer, int rowOrColIndex) {
    node[layer][rowOrColIndex][treeNodeIndex].lzTag = 0;
    if(lowerBound == upperBound) {
        if(layer & 1)
            node[layer][rowOrColIndex][treeNodeIndex].minValue = graph.grids[lowerBound][rowOrColIndex][layer].get_remaining();
        else
            node[layer][rowOrColIndex][treeNodeIndex].minValue = graph.grids[rowOrColIndex][lowerBound][layer].get_remaining();
        return;
    }
    int mid = (lowerBound + upperBound) >> 1;
    build(treeNodeIndex<<1, lowerBound, mid, layer, rowOrColIndex);
    build(treeNodeIndex<<1|1, mid+1, upperBound, layer, rowOrColIndex);
    pushup(treeNodeIndex, layer, rowOrColIndex);
}
//Minimum value in the interval
void SegmentTree::pushup(int treeNodeIndex, int layer, int rowOrColIndex) {
    node[layer][rowOrColIndex][treeNodeIndex].minValue = min(node[layer][rowOrColIndex][treeNodeIndex<<1].minValue, node[layer][rowOrColIndex][treeNodeIndex<<1|1].minValue);
}
//Pushdown lazytag
void SegmentTree::pushdown(int treeNodeIndex, int layer, int rowOrColIndex) {
    if (!node[layer][rowOrColIndex][treeNodeIndex].lzTag) {
        return;
    }
    node[layer][rowOrColIndex][treeNodeIndex<<1].minValue += node[layer][rowOrColIndex][treeNodeIndex].lzTag;
    node[layer][rowOrColIndex][treeNodeIndex<<1|1].minValue += node[layer][rowOrColIndex][treeNodeIndex].lzTag;

    node[layer][rowOrColIndex][treeNodeIndex<<1].lzTag += node[layer][rowOrColIndex][treeNodeIndex].lzTag;
    node[layer][rowOrColIndex][treeNodeIndex<<1|1].lzTag += node[layer][rowOrColIndex][treeNodeIndex].lzTag;
    node[layer][rowOrColIndex][treeNodeIndex].lzTag = 0;
}
//Trigger query function
int SegmentTree::get_remaining_supply(int startIndex, int endIndex, int layer, int rowOrColIndex) {
    if(layer & 1)
        return query(1, 0, graph.row-1, startIndex, endIndex, layer, rowOrColIndex);
    else
        return query(1, 0, graph.column-1, startIndex, endIndex, layer, rowOrColIndex);
}
//Trigger update event
void SegmentTree::update_remaining_supply(int startIndex, int endIndex, int layer, int rowOrColIndex, int delta) {
    if(layer & 1)
        update(1, 0, graph.row-1, startIndex, endIndex, layer, rowOrColIndex, delta);
    else
        update(1, 0, graph.column-1, startIndex, endIndex, layer, rowOrColIndex, delta);
}
//Get the lowest node containing this segment
int SegmentTree::query(int treeNodeIndex, int lowerBound, int upperBound, int startIndex, int endIndex, int layer, int rowOrColIndex) {
    if(startIndex <= lowerBound && endIndex >= upperBound) {
        return node[layer][rowOrColIndex][treeNodeIndex].minValue;
    }
    int mid = (lowerBound + upperBound) >> 1;
    pushdown(treeNodeIndex, layer, rowOrColIndex);
    if(startIndex <= mid) {
        return query(treeNodeIndex<<1, lowerBound, mid, startIndex, endIndex, layer, rowOrColIndex);
    }
    else if(endIndex > mid) {
        return query(treeNodeIndex<<1|1, mid+1, upperBound, startIndex, endIndex, layer, rowOrColIndex);
    }
    else {
        return min(query(treeNodeIndex<<1, lowerBound, mid, startIndex, endIndex, layer, rowOrColIndex), query(treeNodeIndex<<1|1, mid+1, upperBound, startIndex, endIndex, layer, rowOrColIndex));
    }
}
//update segment tree
void SegmentTree::update(int treeNodeIndex, int lowerBound, int upperBound, int startIndex, int endIndex, int layer, int rowOrColIndex, int delta) {
    if(startIndex <= lowerBound && endIndex >= upperBound) {
        node[layer][rowOrColIndex][treeNodeIndex].minValue += delta;
        node[layer][rowOrColIndex][treeNodeIndex].lzTag += delta;
        return;
    }
    int mid = (lowerBound + upperBound) >> 1;
    pushdown(layer, rowOrColIndex, treeNodeIndex);
    if(startIndex <= mid) {
        update(treeNodeIndex<<1, lowerBound, mid, startIndex, endIndex, layer, rowOrColIndex, delta);
    }
    if(endIndex > mid) {
        update(treeNodeIndex<<1|1, mid+1, upperBound, startIndex, endIndex, layer, rowOrColIndex, delta);
    }
    pushup(layer, rowOrColIndex, treeNodeIndex);
}
