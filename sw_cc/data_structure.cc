#include<bits/stdc++.h>
#include "data_structure.h"
#include "segment_tree.h"
using namespace std;

void Net::convert_seg_to_2pin(vector<vector<vector<bool>>>& passingMap, 
        std::vector<Cell>& cellInstances, 
        std::vector<MasterCell>& masterCells
        ) {
    // add segment in passingMap
    for(const auto& segment : routingSegments) {
        if(segment.first.x != segment.second.x) {
            int start = min(segment.first.x, segment.second.x);
            int end = max(segment.first.x, segment.second.x);
            for(int n=start; n<=end; n++) {
                passingMap[n][segment.first.y][segment.first.z];
            }
        }
        else if(segment.first.y != segment.second.y) {
            int start = min(segment.first.y, segment.second.y);
            int end = max(segment.first.y, segment.second.y);
            for(int n=start; n<=end; n++) {
                passingMap[segment.first.x][n][segment.first.z];
            }
        } else if(segment.first.z != segment.second.z) {
            int start = min(segment.first.z, segment.second.z);
            int end = max(segment.first.z, segment.second.z);
            for(int n=start; n<=end; n++) {
                passingMap[segment.first.x][segment.first.y][n];
            }
        }
    }
    // construct pin set 
    unordered_set <Point> pin_map;
    for(const auto& pin:pins) {
        int cell_idx = pin.first;
        int pin_idx = pin.second;
        auto& cell = cellInstances[cell_idx];
        int layer = masterCells[cell.masterCell].pins[pin_idx].layer;
        pin_map.emplace(cell.x,cell.y,layer);    
    }
    // construct 2pin net
    Point start = *pin_map.begin();
}

void Net::traverse_passing_map(vector<vector<vector<bool>>>& passingMap, 
        unordered_set <Point>& pin_map, Point start
        ) {
    
}

int Net::return_node_degree(vector<vector<vector<bool>>>& passingMap, Point p) {
    int x_bound = passingMap.size();
    int y_bound = passingMap[0].size();
    int z_bound = passingMap[0][0].size();
    int count = 0;
    int dir[6] = {1,-1,1,-1,1,-1};
    for(int n=0; n<6; n++) {
        if(n==0 && p.x+dir[n] < x_bound && passingMap[p.x+dir[n]][p.y][p.z])
            count++;
        else if(n==1 && p.x+dir[n] >= 0 && passingMap[p.x+dir[n]][p.y][p.z])
            count++;
        else if(n==2 && p.y+dir[n] < y_bound && passingMap[p.x][p.y+dir[n]][p.z])
            count++;
        else if(n==3 && p.y+dir[n] >= 0 && passingMap[p.x][p.y+dir[n]][p.z])
            count++;
        else if(n==4 && p.z+dir[n] < z_bound && passingMap[p.x][p.y][p.z+dir[n]])
            count++;
        else if(n==5 && p.z+dir[n] >= 0 && passingMap[p.x][p.y][p.z+dir[n]])
            count++;
    }
    return count;
}

RoutingGraph::RoutingGraph() {segmentTree = new SegmentTree(*this);}

RoutingGraph::~RoutingGraph() {delete segmentTree;}

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
    for(auto net_index : cellInstances[cellIndex].connectedNets) {
        del_net_from_graph(net_index);
    }
}

void RoutingGraph::add_net_demand_into_graph(int x, int y, int z, int netIndex) {
    auto hint = grids[x][y][z].passingNets.find(netIndex);
    if(hint == grids[x][y][z].passingNets.end()) {
        grids[x][y][z].demand++;
        grids[x][y][z].passingNets.insert(hint, netIndex);
    }
}
void RoutingGraph::del_net_from_graph(int netIndex) {
    auto& net = nets[netIndex];
    while(net.routingSegments.size() > 0) {
        auto segment = net.routingSegments.back();
        del_seg_demand(segment, netIndex);
        net.routingSegments.pop_back();
    }
}

void RoutingGraph::del_seg_demand(std::pair<Point,Point> segment, int netIndex) {
    int startRow = segment.first.x, startColumn = segment.first.y, startLayer = segment.first.z;
    int endRow = segment.second.x, endColumn = segment.second.y, endLayer = segment.second.z;
    // handle vertical segment
    if(startRow != endRow) {
        for(int j = startRow; j <= endRow; j++) {
            del_seg_demand_from_graph(j, startColumn, startLayer, netIndex);
        }
    }
    // handle horizontal segment
    else if(startColumn != endColumn) {
        for(int j = startColumn; j <= endColumn; j++) {
            del_seg_demand_from_graph(startRow, j, startLayer, netIndex);
        }
    }
    // handle via segment
    else if(startLayer != endLayer) {
        for(int j = startLayer; j <= endLayer; j++) {
            del_seg_demand_from_graph(startRow, startColumn, j, netIndex);
        }
    }
}

void RoutingGraph::del_seg_demand_from_graph(int x, int y, int z, int netIndex) {
    auto hint = grids[x][y][z].passingNets.find(netIndex);
    if(hint != grids[x][y][z].passingNets.end()) {
        grids[x][y][z].demand--;
        grids[x][y][z].passingNets.erase(hint);
    }
}

void RoutingGraph::construct_2pin_nets() {
    // initial pass_map
    vector<vector<vector<bool>>> passingMap;
    passingMap.resize(row);
    for(int r=0; r<row; r++) {
        passingMap[r].resize(column);
        for(int c=0; c<column; c++) {
            passingMap[r][c].resize(layer,0);
        }
    }
    // mark segment passing
    for(const auto& net : nets) {
        
    }

}
