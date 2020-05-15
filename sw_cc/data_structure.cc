#include<bits/stdc++.h>
#include "data_structure.h"

using namespace std;

void RoutingGraph::add_cell_demand_into_graph(int x, int y, int MCtype) {
    //Add blockage demand
    MasterCell& masterCell = masterCells[MCtype];
    for(auto it = masterCell.blockages.begin(); it != masterCell.blockages.end(); ++it) {
        grids[x][y][it->first].demand += it->second;
    }
    //Add sameGGrid rule
    for(auto it = sameGGrid[MCtype].begin(); it != sameGGrid[MCtype].end(); ++it) {
        int MCtype2 = it->first;
        int originalPairCnt = min(cellCount[x][y][MCtype], cellCount[x][y][MCtype2]);
        int afterPairCnt = min(cellCount[x][y][MCtype] + 1, cellCount[x][y][MCtype2]);
        for(auto _it = it->second.begin(); _it != it->second.end(); ++_it) {
            grids[x][y][_it->first].demand += _it->second * (afterPairCnt - originalPairCnt);
        }	
    }
    //Add adjHGGrid rule
    for(auto it = adjHGGrid[MCtype].begin(); it != adjHGGrid[MCtype].end(); ++it) {
        int MCtype2 = it->first, preOriginalPairCnt = 0, preAfterPairCnt = 0, nxtOriginalPairCnt = 0, nxtAfterPairCnt = 0;
        //Update left
        if(y > 0) {
            preOriginalPairCnt = min(cellCount[x][y-1][MCtype2], cellCount[x][y][MCtype]);
            preAfterPairCnt = min(cellCount[x][y-1][MCtype2], cellCount[x][y][MCtype] + 1);
            for(auto _it = it->second.begin(); _it != it->second.end(); ++_it) {
                grids[x][y-1][_it->first].demand += _it->second * (preAfterPairCnt - preOriginalPairCnt);
            }
        }
        //Update right
        if(y < (column - 1)) {
            nxtOriginalPairCnt = min(cellCount[x][y+1][MCtype2], cellCount[x][y][MCtype]);
            nxtAfterPairCnt = min(cellCount[x][y+1][MCtype2], cellCount[x][y][MCtype] + 1);
            for(auto _it = it->second.begin(); _it != it->second.end(); ++_it) {
                grids[x][y+1][_it->first].demand += _it->second * (nxtAfterPairCnt - nxtOriginalPairCnt);
            }
        }
        //Update current
        for(auto _it = it->second.begin(); _it != it->second.end(); ++_it) {
            if((preAfterPairCnt - preOriginalPairCnt + nxtAfterPairCnt - nxtOriginalPairCnt)>0) cout << x << y <<endl;
            grids[x][y][_it->first].demand += _it->second * (preAfterPairCnt - preOriginalPairCnt + nxtAfterPairCnt - nxtOriginalPairCnt);
        }
    }
    cellCount[x][y][MCtype]++;
}

void RoutingGraph::del_cell_demand_from_graph(int x, int y, int MCtype) {
    //Delete blockage demand
    MasterCell& masterCell = masterCells[MCtype];
    for(auto it = masterCell.blockages.begin(); it != masterCell.blockages.end(); ++it) {
        grids[x][y][it->first].demand -= it->second;
    }
    //Delete sameGGrid rule
    for(auto it = sameGGrid[MCtype].begin(); it != sameGGrid[MCtype].end(); ++it) {
        int MCtype2 = it->first;
        int originalPairCnt = min(cellCount[x][y][MCtype], cellCount[x][y][MCtype2]);
        int afterPairCnt = min(cellCount[x][y][MCtype] - 1, cellCount[x][y][MCtype2]);
        for(auto _it = it->second.begin(); _it != it->second.end(); ++_it) {
            grids[x][y][_it->first].demand -= _it->second * (originalPairCnt - afterPairCnt);
        }
    }
    //Delete adjHGGrid rule
    for(auto it = adjHGGrid[MCtype].begin(); it != adjHGGrid[MCtype].end(); ++it) {
        int MCtype2 = it->first, preOriginalPairCnt = 0, preAfterPairCnt = 0, nxtOriginalPairCnt = 0, nxtAfterPairCnt = 0;
        //Delete left
        if(y > 0) {
            preOriginalPairCnt = min(cellCount[x][y-1][MCtype2], cellCount[x][y][MCtype]);
            preAfterPairCnt = min(cellCount[x][y-1][MCtype2], cellCount[x][y][MCtype] - 1);
            for(auto _it = it->second.begin(); _it != it->second.end(); ++_it) {
                grids[x][y-1][_it->first].demand -= _it->second * (preOriginalPairCnt - preAfterPairCnt);
            }
        }
        //Delete right
        if(y < (column - 1)) {
            nxtOriginalPairCnt = min(cellCount[x][y+1][MCtype2], cellCount[x][y][MCtype]);
            nxtAfterPairCnt = min(cellCount[x][y+1][MCtype2], cellCount[x][y][MCtype] - 1);
            for(auto _it = it->second.begin(); _it != it->second.end(); ++_it) {
                grids[x][y+1][_it->first].demand -= _it->second * (nxtOriginalPairCnt - nxtAfterPairCnt);
            }
        }
        //Delete current
        for(auto _it = it->second.begin(); _it != it->second.end(); ++_it) {
            grids[x][y][_it->first].demand -= _it->second * (preOriginalPairCnt  - preAfterPairCnt + nxtOriginalPairCnt - nxtAfterPairCnt);
        }
    }
    cellCount[x][y][MCtype]--;
}

void RoutingGraph::add_cell(int x, int y, int cellIndex) {
    //Add into placement
    placement[x][y].insert(cellIndex);
    //Add demand into graph
    add_cell_demand_into_graph(x, y, cellInstances[cellIndex].masterCell);
}

void RoutingGraph::del_cell(int cellIndex) {
    Cell& cell = cellInstances[cellIndex];
    int x = cell.x;
    int y = cell.y;
    //Del from placement
    placement[x][y].erase(cellIndex);
    //Del demand from graph
    del_cell_demand_from_graph(x, y, cellInstances[cellIndex].masterCell);
}

void RoutingGraph::add_net_demand_into_graph(int x, int y, int z, int netIndex) {
    auto hint = grids[x][y][z].passingNets.find(netIndex);
    if(hint == grids[x][y][z].passingNets.end()) {
        grids[x][y][z].demand++;
        grids[x][y][z].passingNets.insert(hint, netIndex);
    }
}
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
//Trigger query function
int SegmentTree::get_remaining_supply(int startIndex, int endIndex, int layer, int rowOrColIndex) {
    if(layer & 1) 
        return query(1, 0, graph.row-1, startIndex, endIndex, layer, rowOrColIndex);
    else 
        return query(1, 0, graph.column-1, startIndex, endIndex, layer, rowOrColIndex);
}
//Get the lowest node containing this segment
int SegmentTree::query(int treeNodeIndex, int lowerBound, int upperBound, int startIndex, int endIndex, int layer, int rowOrColIndex) {
    if(startIndex <= lowerBound && endIndex >= upperBound) {
        return node[layer][rowOrColIndex][treeNodeIndex].minValue;
    }
    int mid = (lowerBound + upperBound) >> 1;
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

