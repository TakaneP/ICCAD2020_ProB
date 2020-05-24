#include<bits/stdc++.h>
#include "data_structure.h"
#include "segment_tree.h"
using namespace std;

void Net::convert_seg_to_2pin(vector<vector<vector<DegreeNode>>>& degreeMap, 
        std::vector<Cell>& cellInstances, 
        std::vector<MasterCell>& masterCells
        ) {
    // construct pin set 
    unordered_set <Point,MyHashFunction> pin_map;
    unordered_set <Point,MyHashFunction> steiner_map;
    for(const auto& pin:pins) {
        int cell_idx = pin.first;
        int pin_idx = pin.second;
        auto& cell = cellInstances[cell_idx];
        int layer = masterCells[cell.masterCell].pins[pin_idx].layer;
        Point p(cell.x,cell.y,layer);
        if(pin_map.find(p) != pin_map.end()) {
            // pin in same GCell
            if (minRoutingLayer >= 0) {
                steiner_map.emplace(cell.x,cell.y,minRoutingLayer);    
            }
            TwoPinNet same_cell_net;
            same_cell_net.n1.p = p;
            same_cell_net.n1.type = 1;
            same_cell_net.n2.p = p;
            same_cell_net.n2.type = 1;
            same_cell_net.paths.push_back(std::pair<Point,Point>(p,p));
            routingTree.push_back(same_cell_net);
        }
        else
            pin_map.emplace(cell.x,cell.y,layer);    
    }
    // add segment in passingMap, and construct steiner_map
    for(const auto& segment : routingSegments) {
        if(segment.first.x != segment.second.x) {
            int start = min(segment.first.x, segment.second.x);
            int end = max(segment.first.x, segment.second.x);
            for(int n=start; n<=end; n++) {
                if(n == start)
                    degreeMap[n][segment.first.y][segment.first.z].right = 1;
                else if(n == end)
                    degreeMap[n][segment.first.y][segment.first.z].left = 1;
                else {
                    degreeMap[n][segment.first.y][segment.first.z].right = 1;
                    degreeMap[n][segment.first.y][segment.first.z].left = 1;
                }
                if(degreeMap[n][segment.first.y][segment.first.z].return_degree() >= 3) 
                    steiner_map.emplace(n,segment.first.y,segment.first.z);
            }       
        }
        else if(segment.first.y != segment.second.y) {
            int start = min(segment.first.y, segment.second.y);
            int end = max(segment.first.y, segment.second.y);
            for(int n=start; n<=end; n++) {
                if(n == start)
                    degreeMap[segment.first.x][n][segment.first.z].up = 1;
                else if(n == end)
                    degreeMap[segment.first.x][n][segment.first.z].down = 1;
                else {
                    degreeMap[segment.first.x][n][segment.first.z].up = 1;
                    degreeMap[segment.first.x][n][segment.first.z].down = 1;
                }
                if(degreeMap[segment.first.x][n][segment.first.z].return_degree() >= 3)
                    steiner_map.emplace(segment.first.x,n,segment.first.z);
            }
        } else if(segment.first.z != segment.second.z) {
            int start = min(segment.first.z, segment.second.z);
            int end = max(segment.first.z, segment.second.z);
            for(int n=start; n<=end; n++) {
                if(n == start)
                    degreeMap[segment.first.x][segment.first.y][n].top = 1;
                else if(n == end)
                    degreeMap[segment.first.x][segment.first.y][n].bottom = 1;
                else {
                    degreeMap[segment.first.x][segment.first.y][n].top = 1;
                    degreeMap[segment.first.x][segment.first.y][n].bottom = 1;
                }
                if(degreeMap[segment.first.x][segment.first.y][n].return_degree() >= 3) 
                    steiner_map.emplace(segment.first.x,segment.first.y,n);
            }
        }
    }
    // construct 2pin net
    Point start = *pin_map.begin();
    traverse_passing_map(degreeMap, pin_map, steiner_map, start);
    print_two_pins();
}

void Net::traverse_passing_map(vector<vector<vector<DegreeNode>>>& degreeMap, 
        unordered_set <Point,MyHashFunction>& pin_map, 
        unordered_set <Point,MyHashFunction>& steiner_map, 
        Point start_p
        ) {
    bool is_pin = (pin_map.find(start_p) != pin_map.end());
    while(true) {  // for 6 dir     
        Point dir = return_next_dir(degreeMap, start_p);
        if(dir == Point(-1,-1,-1))
            break;
        // -dir because delete outgoing edge
        decrese_degree_endpoint(degreeMap, start_p, Point(-dir.x, -dir.y, -dir.z));
        Point now_p = start_p;
        TwoPinNet segment;
        segment.n1.p = start_p;
        segment.n1.type = (is_pin) ? 1 : 0;
        std::pair<Point,Point> path;
        path.first = start_p;
        while(true) { // for traversing
            now_p = now_p+dir;
            bool is_pin = pin_map.find(now_p) != pin_map.end();
            bool is_steiner = steiner_map.find(now_p) != steiner_map.end();
            if(is_pin || is_steiner) {
                // create 2pin
                decrese_degree_endpoint(degreeMap, now_p, dir);
                path.second = now_p;
                segment.paths.push_back(path);
                segment.n2.p = now_p;
                segment.n2.type = (is_pin) ? 1 : 0;
                this->routingTree.push_back(segment);              
                traverse_passing_map(degreeMap, pin_map, steiner_map, now_p);
                break;
            }
            else {
                // keep propogate
                Point next_p = now_p + dir;            
                if(check_map_legal(degreeMap, next_p) && check_map_dir(degreeMap, now_p, dir)) {
                    decrese_degree_middle_p(degreeMap, now_p, dir);
                }
                else {
                    // create segment path
                    decrese_degree_endpoint(degreeMap, now_p, dir);
                    path.second = now_p;
                    segment.paths.push_back(path);
                    path.first = now_p;
                    dir = return_next_dir(degreeMap, now_p);
                    if(dir == Point(-1,-1,-1)) {
                        // redundant net
                        segment.n2.p = now_p;
                        segment.n2.type = -1;
                        this->routingTree.push_back(segment); 
                        break;
                    } else {
                        decrese_degree_endpoint(degreeMap, now_p, Point(-dir.x, -dir.y, -dir.z));
                    }                   
                }
            }
        }
    }
}

Point Net::return_next_dir(vector<vector<vector<DegreeNode>>>& degreeMap, Point now_p) {
    if(degreeMap[now_p.x][now_p.y][now_p.z].right)
        return Point(1,0,0);
    else if(degreeMap[now_p.x][now_p.y][now_p.z].left)
        return Point(-1,0,0);
    else if(degreeMap[now_p.x][now_p.y][now_p.z].up)
        return Point(0,1,0);
    else if(degreeMap[now_p.x][now_p.y][now_p.z].down)
        return Point(0,-1,0);
    else if(degreeMap[now_p.x][now_p.y][now_p.z].top)
        return Point(0,0,1);
    else if(degreeMap[now_p.x][now_p.y][now_p.z].bottom)
        return Point(0,0,-1);
    else
        return Point(-1,-1,-1);
}

bool Net::check_map_legal(vector<vector<vector<DegreeNode>>>& degreeMap, Point now_p) {
    int x_bound = degreeMap.size();
    int y_bound = degreeMap[0].size();
    int z_bound = degreeMap[0][0].size();
    return (now_p.x<x_bound && now_p.x>=0 && 
            now_p.y<y_bound && now_p.y>=0 && 
            now_p.z<z_bound && now_p.z>=0);
}

bool Net::check_map_dir(vector<vector<vector<DegreeNode>>>& degreeMap, Point now_p, Point dir) {
    if(dir.x==1)
        return degreeMap[now_p.x][now_p.y][now_p.z].right; 
    if(dir.x==-1)
        return degreeMap[now_p.x][now_p.y][now_p.z].left;
    if(dir.y==1)
        return degreeMap[now_p.x][now_p.y][now_p.z].up;
    if(dir.y==-1)       
        return degreeMap[now_p.x][now_p.y][now_p.z].down;
    if(dir.z==1)
        return degreeMap[now_p.x][now_p.y][now_p.z].top;
    if(dir.z==-1)
        return degreeMap[now_p.x][now_p.y][now_p.z].bottom;       
}

void Net::decrese_degree_endpoint(std::vector<std::vector<std::vector<DegreeNode>>>& degreeMap, 
        Point now_p, Point dir) {
    if(dir.x==1)
        degreeMap[now_p.x][now_p.y][now_p.z].left = 0;
    if(dir.x==-1)
        degreeMap[now_p.x][now_p.y][now_p.z].right = 0;
    if(dir.y==1)
        degreeMap[now_p.x][now_p.y][now_p.z].down = 0;
    if(dir.y==-1)
        degreeMap[now_p.x][now_p.y][now_p.z].up = 0;
    if(dir.z==1)
        degreeMap[now_p.x][now_p.y][now_p.z].bottom = 0;
    if(dir.z==-1)
        degreeMap[now_p.x][now_p.y][now_p.z].top = 0;
}

void Net::decrese_degree_middle_p(std::vector<std::vector<std::vector<DegreeNode>>>& degreeMap, 
        Point now_p, Point dir) {
    if(dir.x!=0) {
        degreeMap[now_p.x][now_p.y][now_p.z].left = 0;
        degreeMap[now_p.x][now_p.y][now_p.z].right = 0;
    }
    if(dir.y!=0) {
        degreeMap[now_p.x][now_p.y][now_p.z].up = 0;
        degreeMap[now_p.x][now_p.y][now_p.z].down = 0;
    }
    if(dir.z!=0) {
        degreeMap[now_p.x][now_p.y][now_p.z].top = 0;
        degreeMap[now_p.x][now_p.y][now_p.z].bottom = 0;
    }
}

void Net::print_two_pins() {
    for(auto twopin : this->routingTree) {
        cout << endl << twopin.n1.p << " " << twopin.n1.type << " to " << 
        twopin.n2.p << " " << twopin.n2.type << "\n";
        for(auto path : twopin.paths) {
            cout << "path: " << path.first << " " << path.second << endl;
        }
    }  
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
    vector<vector<vector<DegreeNode>>> degreeMap;
    degreeMap.resize(row);
    for(int r=0; r<row; r++) {
        degreeMap[r].resize(column);
        for(int c=0; c<column; c++) {
            degreeMap[r][c].resize(layer);
        }
    }
    // mark segment passing
    int a=0;
    for(auto& net : nets) {
        cout << "\nNew net: " << a++ << "\n";
        net.convert_seg_to_2pin(degreeMap, cellInstances, masterCells);
    }

}
