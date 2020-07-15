#include<bits/stdc++.h>
#include "data_structure.h"
#include "segment_tree.h"

using namespace std;

bool operator>(const TwoPinNet& n1, const TwoPinNet& n2) {
    return n1.wire_length>n2.wire_length;
}

size_t MyHashFunction::operator()(const Point& p) const {
    return (p.x+p.y*2000+p.z*2000*2000);
}

void Net::convert_seg_to_2pin(vector<vector<vector<DegreeNode>>>& degreeMap, 
        std::vector<Cell>& cellInstances, 
        std::vector<MasterCell>& masterCells
        ) {
    // construct pin set 
    unordered_set <Point,MyHashFunction> pin_map;
    unordered_set <Point,MyHashFunction> steiner_map;
    int row = degreeMap.size();
    int column = (degreeMap.empty())? 0:degreeMap[0].size();
    vector<bool> checkRedundant(row*column, false);
    map<tuple<int,int,int>, vector<pair<int,int>>> localNets;
    for(const auto& pin:pins) {
        int cell_idx = pin.first;
        int pin_idx = pin.second;
        auto& cell = cellInstances[cell_idx];
        int layer = cell.pins[pin_idx].layer;
        Point p(cell.x,cell.y,layer);
        localNets[{cell.x, cell.y, layer}].push_back({cell_idx, pin_idx});
        pin_map.emplace(cell.x,cell.y,layer);
        /*
        if(pin_map.find(p) != pin_map.end()) {
            // pin in same GCell
            if (minRoutingLayer > layer && !checkRedundant[cell.x*row + cell.y]) {
                steiner_map.emplace(cell.x,cell.y,minRoutingLayer);
                checkRedundant[cell.x*row + cell.y] = true;   
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
            pin_map.emplace(cell.x,cell.y,layer);*/
    }
    // add segment in passingMap, and construct steiner_map
    set_passing_map(degreeMap, cellInstances, masterCells, pin_map, steiner_map, 1);
    // construct 2pin net
    Point start = *pin_map.begin();
    std::vector<TwoPinNet> routingTree;
    traverse_passing_map(degreeMap, pin_map, steiner_map, start, localNets, routingTree);
    // reset passingMap
    set_passing_map(degreeMap, cellInstances, masterCells, pin_map, steiner_map, 0);
    //print_two_pins();
    construct_branch_nodes(routingTree);
    remove_dangling_wire();
    //print_two_pins();
    remove_branch_cycle();
}

void Net::set_passing_map(vector<vector<vector<DegreeNode>>>& degreeMap, std::vector<Cell>& cellInstances, std::vector<MasterCell>& masterCells, 
    unordered_set <Point,MyHashFunction>& pin_map, unordered_set <Point,MyHashFunction>& steiner_map, int value) {
    // add segment in passingMap, and construct steiner_map
    for(const auto& segment : routingSegments) {
        if(segment.first.x != segment.second.x) {
            int start = min(segment.first.x, segment.second.x);
            int end = max(segment.first.x, segment.second.x);
            for(int n=start; n<=end; n++) {
                if(n == start)
                    degreeMap[n][segment.first.y][segment.first.z].right = value;
                else if(n == end)
                    degreeMap[n][segment.first.y][segment.first.z].left = value;
                else {
                    degreeMap[n][segment.first.y][segment.first.z].right = value;
                    degreeMap[n][segment.first.y][segment.first.z].left = value;
                }
                if(value && degreeMap[n][segment.first.y][segment.first.z].return_degree() >= 3) 
                    steiner_map.emplace(n,segment.first.y,segment.first.z);
            }       
        }
        else if(segment.first.y != segment.second.y) {
            int start = min(segment.first.y, segment.second.y);
            int end = max(segment.first.y, segment.second.y);
            for(int n=start; n<=end; n++) {
                if(n == start)
                    degreeMap[segment.first.x][n][segment.first.z].up = value;
                else if(n == end)
                    degreeMap[segment.first.x][n][segment.first.z].down = value;
                else {
                    degreeMap[segment.first.x][n][segment.first.z].up = value;
                    degreeMap[segment.first.x][n][segment.first.z].down = value;
                }
                if(value && degreeMap[segment.first.x][n][segment.first.z].return_degree() >= 3)
                    steiner_map.emplace(segment.first.x,n,segment.first.z);
            }
        } else if(segment.first.z != segment.second.z) {
            int start = min(segment.first.z, segment.second.z);
            int end = max(segment.first.z, segment.second.z);
            for(int n=start; n<=end; n++) {
                if(n == start)
                    degreeMap[segment.first.x][segment.first.y][n].top = value;
                else if(n == end)
                    degreeMap[segment.first.x][segment.first.y][n].bottom = value;
                else {
                    degreeMap[segment.first.x][segment.first.y][n].top = value;
                    degreeMap[segment.first.x][segment.first.y][n].bottom = value;
                }
                if(value && degreeMap[segment.first.x][segment.first.y][n].return_degree() >= 3) 
                    steiner_map.emplace(segment.first.x,segment.first.y,n);
            }
        }
    }
}

void Net::traverse_passing_map(vector<vector<vector<DegreeNode>>>& degreeMap, 
        unordered_set <Point,MyHashFunction>& pin_map, 
        unordered_set <Point,MyHashFunction>& steiner_map, 
        Point start_p,
        map<tuple<int,int,int>, vector<pair<int,int>>>& localNets,
        std::vector<TwoPinNet>& routingTree
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
        if(localNets[{now_p.x, now_p.y, now_p.z}].size() > 1) {
            segment.n1.type = 2;
            segment.n1.mergedLocalPins = localNets[{now_p.x, now_p.y, now_p.z}];
        }
        else
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
                if(localNets[{now_p.x, now_p.y, now_p.z}].size() > 1) {
                    segment.n2.type = 2;
                    segment.n2.mergedLocalPins = localNets[{now_p.x, now_p.y, now_p.z}];
                }
                else
                    segment.n2.type = (is_pin) ? 1 : 0;
                if(!(segment.n1.p == segment.n2.p))
                    routingTree.push_back(segment);              
                traverse_passing_map(degreeMap, pin_map, steiner_map, now_p, localNets, routingTree);
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
                        routingTree.push_back(segment); 
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

void Net::print_two_pins(std::vector<TwoPinNet>& routingTree) {
    for(auto twopin : routingTree) {
        cout << endl << twopin.n1.p << " " << twopin.n1.type << " to " << 
        twopin.n2.p << " " << twopin.n2.type << "\n";
        for(auto path : twopin.paths) {
            cout << "path: " << path.first << " " << path.second << endl;
        }
        if(twopin.n1.type == 2) {
            cout << "local net Node1: \n";
            for(auto& pin : twopin.n1.mergedLocalPins) cout << "Cell " << pin.first+1 << " Pin " << pin.second+1  << " ";
            cout << "\n";
        }
        if(twopin.n2.type == 2) {
            cout << "local net Node2: \n";
            for(auto& pin : twopin.n2.mergedLocalPins) cout << "Cell " << pin.first+1 << " Pin " << pin.second+1 << " ";
            cout << "\n";
        }
        twopin.update_wire_length();
        cout << "wl: " << twopin.wire_length << endl;
    }  
}

void Net::construct_branch_nodes(std::vector<TwoPinNet>& routingTree) {
    // construct neighbor relation
    for(auto& twopin : routingTree) {
        Point p1 = twopin.n1.p;
        Point p2 = twopin.n2.p;
        auto p1_ptr = branch_nodes.find(p1);
        auto p2_ptr = branch_nodes.find(p2);
        if(p1_ptr == branch_nodes.end()) {
            branch_nodes.emplace(p1, TreeNode(twopin.n1));
            p1_ptr = branch_nodes.find(p1);      
        }
        if(p2_ptr == branch_nodes.end()) {
            branch_nodes.emplace(p2, TreeNode(twopin.n2));
            p2_ptr = branch_nodes.find(p2);
        }
        p1_ptr->second.neighbors.push_back(make_pair(p2,twopin));
        p2_ptr->second.neighbors.push_back(make_pair(p1,twopin));
    }  
}

void Net::remove_dangling_wire() {
    queue<Point> todo_points;
    for(auto& treenode : this->branch_nodes) {
        todo_points.push(treenode.first);
    }
    while(todo_points.size() > 0) {
        Point tree_p = todo_points.front();
        todo_points.pop();
        if (branch_nodes.find(tree_p) == branch_nodes.end())
            continue;
        auto& treenode = branch_nodes[tree_p];
        // find dangling endpoint
        if(treenode.node.type == -1) {
            Point neighbor = treenode.neighbors[0].first;
            TreeNode& effect_node = branch_nodes[neighbor];
            todo_points.push(neighbor);
            for(int n=0; n<effect_node.neighbors.size(); n++) {           
                if(effect_node.neighbors[n].first == tree_p) {
                    // remove dangling endpoint's neighbor to dangling
                    effect_node.neighbors.erase(effect_node.neighbors.begin()+n);               
                    if(effect_node.node.type == 0 && effect_node.neighbors.size()<=1) {
                        effect_node.node.type = -1;
                    }
                    break;
                }
            }
            // remove dangling twopin-net
            branch_nodes.erase(tree_p);
        }
    }
}

void Net::remove_branch_cycle() {
    std::priority_queue<TwoPinNet, vector<TwoPinNet>, greater<TwoPinNet>> frontier_edges;
    std::queue<TwoPinNet> skip_edges;
    std::unordered_set <Point,MyHashFunction> visited_nodes;
    push_edge_in_queue(frontier_edges);
    // construct MST
    while(!frontier_edges.empty()) {
        auto edge = frontier_edges.top();
        frontier_edges.pop();
        //cout << "Two_pin: " << edge.n1.p << " " << edge.n1.type << " to " << 
        //edge.n2.p << " " << edge.n2.type << "\n";
        if(visited_nodes.size() == 0) {
            visited_nodes.insert(edge.n1.p);
            visited_nodes.insert(edge.n2.p);
            continue;
        }
        auto e1_iter = visited_nodes.find(edge.n1.p);
        auto e2_iter = visited_nodes.find(edge.n2.p);
        Point non_tree_node, tree_node;
        if(e1_iter==visited_nodes.end() && e2_iter!=visited_nodes.end()) {
            visited_nodes.insert(edge.n1.p);  
            while(!skip_edges.empty()) {
                frontier_edges.push(skip_edges.front());
                skip_edges.pop();
            }
            continue;
        } 
        if(e2_iter==visited_nodes.end() && e2_iter!=visited_nodes.end()) {
            visited_nodes.insert(edge.n2.p);  
            while(!skip_edges.empty()) {
                frontier_edges.push(skip_edges.front());
                skip_edges.pop();
            }
            continue;
        } 
        if(e1_iter!=visited_nodes.end() && e2_iter!=visited_nodes.end()) {
            //delete this edge in branch_nodes
            auto& e1_treeNode = branch_nodes[edge.n1.p];
            auto& e2_treeNode = branch_nodes[edge.n2.p];
            for(int n=0; n<e1_treeNode.neighbors.size(); n++) {
                if(e1_treeNode.neighbors[n].first == edge.n2.p) {
                    e1_treeNode.neighbors.erase(e1_treeNode.neighbors.begin()+n);
                    break;
                }
            }
            for(int n=0; n<e2_treeNode.neighbors.size(); n++) {
                if(e2_treeNode.neighbors[n].first == edge.n1.p) {
                    e2_treeNode.neighbors.erase(e2_treeNode.neighbors.begin()+n);
                    break;
                }
            }
            continue;
        } else {
            skip_edges.push(edge);
        }
    }
}

void Net::push_edge_in_queue(std::priority_queue<TwoPinNet, vector<TwoPinNet>, greater<TwoPinNet>>& frontier_edges) {
    // push edge in frontier_edges
    std::unordered_set <Point,MyHashFunction> visited_points;
    std::unordered_set <Point,MyHashFunction> traverse_points;
    if(branch_nodes.size() == 0) return;
    Point bfs_p = branch_nodes.begin()->first;
    traverse_points.insert(bfs_p);
    while(!traverse_points.empty()) { 
        bfs_p = *traverse_points.begin();
        traverse_points.erase(bfs_p);
        auto& treenode = branch_nodes[bfs_p];
        // use map to solve multi-edge, and record less wire length edge
        unordered_map<Point,TwoPinNet,MyHashFunction> edges;
        for(auto& neighbor : treenode.neighbors) {                     
            if(visited_points.find(neighbor.first) != visited_points.end()) {
                // visited before
                continue;
            }
            neighbor.second.update_wire_length();
            if(edges.find(neighbor.first) != edges.end()) {
                if(edges[neighbor.first].wire_length > neighbor.second.wire_length)
                    edges[neighbor.first] = neighbor.second;
            }         
        }
        for(auto iter : edges) {
            frontier_edges.push(iter.second);
            traverse_points.insert(iter.first);
        }
        visited_points.insert(bfs_p);
    }
}

RoutingGraph::RoutingGraph(): usedCellMove(0) {segmentTree = new SegmentTree(*this);}
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
    cellInstances[cellIndex].x = x;
    cellInstances[cellIndex].y = y;
    //Add demand into graph
    add_cell_demand_into_graph(x, y, cellInstances[cellIndex].mcType);
}

void RoutingGraph::del_cell(int cellIndex) {
    Cell& cell = cellInstances[cellIndex];
    int x = cell.x;
    int y = cell.y;
    //Del from placement
    placement[x][y].erase(cellIndex);
    //Del demand from graph
    del_cell_demand_from_graph(x, y, cellInstances[cellIndex].mcType);
    for(auto pin : cellInstances[cellIndex].pins) {
        int netIndex = pin.connectedNet;
        if(netIndex != -1) del_net_from_graph(x, y, pin.layer, netIndex);
    }
}

void RoutingGraph::add_net_demand_into_graph(int x, int y, int z, int netIndex) {
    if(grids[x][y][z].passingNets[netIndex]++ == 0) //every pin and segment contribute 1 to passingNet
        grids[x][y][z].demand++;
}
void RoutingGraph::del_net_from_graph(int x, int y, int z, int netIndex) {
    if(--grids[x][y][z].passingNets[netIndex] == 0)
        grids[x][y][z].demand--;
    else {
        auto& net = nets[netIndex];
        while(net.routingSegments.size() > 0) {
            auto segment = net.routingSegments.back();
            del_seg_demand(segment, netIndex);
            net.routingSegments.pop_back();
        }
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
    if(--grids[x][y][z].passingNets[netIndex] == 0)
        grids[x][y][z].demand--;
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

Tree RoutingGraph::RSMT(vector<int> x, vector<int> y) {
    if(x.size() != y.size()) {
        printf("Step RSMT construction: Error, size of x not equal to size of y\n");
        exit(1);
	}
    int d = x.size();
    Tree flutetree;
    int flutewl;
    readLUT();
    int *x_arr = &x[0];
    int *y_arr = &y[0];
    flutetree = flute(d, x_arr, y_arr, ACCURACY);
    printf("FLUTE wirelength = %d\n", flutetree.length);
    flutewl = flute_wl(d, x_arr, y_arr, ACCURACY);
    printf("FLUTE wirelength (without RSMT construction) = %d\n", flutewl);
    return flutetree;
}

bool sortbysec(const pair<Point,int> &a, const pair<Point,int> &b) 
{ 
    return (a.second > b.second); 
} 

void RoutingGraph::move_cells_force() {
    int c=0;
    for(auto& cell : cellInstances) {
        cout << "cell " << c++ << endl;
        if(cell.movable == 0)
            continue;
        // move cell to free space
        vector<int> x_series, y_series;
        for(auto& pin : cell.pins) {
            // find optimal region
            Net& neighbor_net = nets[pin.connectedNet];     
            int min_x=INT_MAX, max_x=0, min_y=INT_MAX, max_y=0;       
            for(auto& net_pin : neighbor_net.pins) {
                int cell_idx = net_pin.first;
                int c_x = cellInstances[cell_idx].x;
                int c_y = cellInstances[cell_idx].y;
                if(c_x == cell.x && c_y == cell.y)
                    continue;  
                min_x = (min_x > c_x) ? c_x : min_x;
                max_x = (max_x < c_x) ? c_x : max_x;
                min_y = (min_y > c_y) ? c_y : min_y;
                max_y = (max_y < c_y) ? c_y : max_y;
            }
            x_series.push_back(min_x);
            x_series.push_back(max_x);
            y_series.push_back(min_y);
            y_series.push_back(max_y);
        }
        // no optimal region
        if(x_series.size()%2 && y_series.size()%2)
            continue;
        sort(x_series.begin(), x_series.end());
        sort(y_series.begin(), y_series.end());
        int opt_x_left = x_series[x_series.size()/2];
        int opt_x_right = x_series[x_series.size()/2+1];
        int opt_y_left = y_series[y_series.size()/2];
        int opt_y_right = y_series[y_series.size()/2+1];
        // no optimal region
        if(opt_x_left == opt_x_right && opt_y_left == opt_y_right)
            continue;
        cout << opt_x_left << " " << opt_y_left << " " << opt_x_right << " " << opt_y_right << endl;
        vector<pair<Point,int>> cells_pos;
        for(int x=opt_x_left; x<=opt_x_right; x++) {
            for(int y=opt_y_left; y<=opt_y_right; y++) {
                int profit = check_cell_cost_in_graph(x, y, cell.mcType);
                cells_pos.emplace_back(Point(x,y,0), profit);
            }
        }
        sort(cells_pos.begin(), cells_pos.end(), sortbysec);
        for(int n=0; n<cells_pos.size(); n++) {
            cout << "(" << cells_pos[n].first.x << ", " << cells_pos[n].first.y << ") cost: " << cells_pos[n].second << endl;
        }
    }
}

int RoutingGraph::check_cell_cost_in_graph(int x, int y, int MCtype) {
    vector<int> layer_remain(layer);
    MasterCell& masterCell = masterCells[MCtype];
    // initial remain
    for(int n=0; n<layer_remain.size(); n++) {
        layer_remain[n] = grids[x][y][n].capacity - grids[x][y][n].demand;
    }
    // add blockage demand
    for(auto it = masterCell.blockages.begin(); it != masterCell.blockages.end(); ++it) {
        layer_remain[it->first] -= it->second;
    }
    // Add sameGGrid rule
    for(auto it = sameGGrid[MCtype].begin(); it != sameGGrid[MCtype].end(); ++it) {
        int MCtype2 = it->first;
        int originalPairCnt = min(cellCount[x][y][MCtype], cellCount[x][y][MCtype2]);
        int afterPairCnt = min(cellCount[x][y][MCtype] + 1, cellCount[x][y][MCtype2]);
        for(auto _it = it->second.begin(); _it != it->second.end(); ++_it) {
            layer_remain[_it->first] -= _it->second * (afterPairCnt - originalPairCnt);
        }	
    }
    // Add adjHGGrid rule
    for(auto it = adjHGGrid[MCtype].begin(); it != adjHGGrid[MCtype].end(); ++it) {
        int MCtype2 = it->first, preOriginalPairCnt = 0, preAfterPairCnt = 0, nxtOriginalPairCnt = 0, nxtAfterPairCnt = 0;
        // check left cell demand
        if(y > 0) {
            preOriginalPairCnt = min(cellCount[x][y-1][MCtype2], cellCount[x][y][MCtype]);
            preAfterPairCnt = min(cellCount[x][y-1][MCtype2], cellCount[x][y][MCtype] + 1);
            for(auto _it = it->second.begin(); _it != it->second.end(); ++_it) {
                int new_demand = grids[x][y-1][_it->first].demand + _it->second * (preAfterPairCnt - preOriginalPairCnt);
                if(grids[x][y-1][_it->first].capacity < new_demand) {
                    return 0;
                }
            }
        }
        // check right cell demand
        if(y < (column - 1)) {
            nxtOriginalPairCnt = min(cellCount[x][y+1][MCtype2], cellCount[x][y][MCtype]);
            nxtAfterPairCnt = min(cellCount[x][y+1][MCtype2], cellCount[x][y][MCtype] + 1);
            for(auto _it = it->second.begin(); _it != it->second.end(); ++_it) {
                int new_demand = grids[x][y+1][_it->first].demand + _it->second * (nxtAfterPairCnt - nxtOriginalPairCnt);
                if(grids[x][y+1][_it->first].capacity < new_demand) {
                    return 0;
                }
            }
        }
        // Update current
        for(auto _it = it->second.begin(); _it != it->second.end(); ++_it) {
            layer_remain[_it->first] -= _it->second * (preAfterPairCnt - preOriginalPairCnt + nxtAfterPairCnt - nxtOriginalPairCnt);
        }
        // calculate profit
        int profit = 0;
        for(int n=0; n<layer_remain.size(); n++) {
            profit += layer_remain[n];
        }
        return profit;
    }
}

void TwoPinNet::update_wire_length() {
    wire_length = 0;
    for(auto& path : paths) {
        wire_length += abs((path.first.x-path.second.x)) + abs((path.first.y-path.second.y)) + abs((path.first.z-path.second.z));
    }
}