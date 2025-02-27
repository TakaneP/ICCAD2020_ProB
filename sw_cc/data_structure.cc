#include<bits/stdc++.h>
#include "data_structure.h"
#include "segment_tree.h"
#include "global_parameter.h"

using namespace std;

void print_neighbors(Net& n) {
    printf("Net %d\n", n.netId);
    for(auto it = n.branch_nodes.begin(); it != n.branch_nodes.end(); ++it) {
        const Point& p = it->first;
        TreeNode& treeNode = it->second;
        printf("Point %d %d %d type %d:\n", p.x, p.y, p.z, treeNode.node.type);
        for(auto _it = treeNode.neighbors.begin(); _it != treeNode.neighbors.end(); _it++) {
            printf("%d %d %d\n", _it->first.x, _it->first.y, _it->first.z);
        }
    }
}

bool operator>(const TwoPinNet& n1, const TwoPinNet& n2) {
    return n1.wire_length>n2.wire_length;
}

bool operator<(const pair<Point,int>& p1, const pair<Point,int>& p2) {
    return p1.second < p2.second;
}

bool operator>(const pair<Point,int>& p1, const pair<Point,int>& p2) {
    return p1.second > p2.second;
}

bool in_range(int target, int a, int b) {
    return target > min(a,b) && target < max(a,b);
}

int distance(const Point& p1, const Point& p2) {
    return abs((p1.x-p2.x)) + abs((p1.y-p2.y)) + abs((p1.z-p2.z));
}

Point norm(const Point& p) {
    return Point( (p.x!=0), (p.y!=0), (p.z!=0) );
}

TwoPinNet two_pin_reverse(TwoPinNet two_pin) {
    swap(two_pin.n1, two_pin.n2);
    reverse(two_pin.paths.begin(), two_pin.paths.end());
    for(auto& path : two_pin.paths)
        swap(path.first, path.second);
    return two_pin;
}

size_t MyHashFunction::operator()(const Point& p) const {
    return (p.x+p.y*2000+p.z*2000*2000);
}

size_t MyHashFunction::operator()(const pair<Point,Point>& p) const {
    return (p.first.x+p.first.y*2000+p.first.z*2000*2000);
}

void Net::convert_seg_to_2pin(vector<vector<vector<DegreeNode>>>& degreeMap, 
        std::vector<Cell>& cellInstances, 
        std::vector<MasterCell>& masterCells,
        vector<vector<vector<Gcell>>>& grids
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
    }
    if(routingSegments.empty()) {
        if(pin_map.empty()) return;
        const Point& localPosition = *(pin_map.begin());
        if(localNets[{localPosition.x, localPosition.y, localPosition.z}].size() == 1) branch_nodes[localPosition].node.type = 1;
        else if(localNets[{localPosition.x, localPosition.y, localPosition.z}].size() > 1) branch_nodes[localPosition].node.type = 2;
        branch_nodes[localPosition].node.mergedLocalPins = localNets[{localPosition.x, localPosition.y, localPosition.z}];
        return;
    }
    // add segment in passingMap, and construct steiner_map
    set_passing_map(degreeMap, cellInstances, masterCells, pin_map, steiner_map, 1);
    // construct 2pin net
    Point start = *pin_map.begin();
    std::vector<TwoPinNet> routingTree;
    traverse_passing_map(degreeMap, pin_map, steiner_map, start, localNets, routingTree, grids);
    // reset passingMap
    set_passing_map(degreeMap, cellInstances, masterCells, pin_map, steiner_map, 0);
    construct_branch_nodes(routingTree);
    remove_dangling_wire(grids);
    //print_two_pins(routingTree);
    remove_branch_cycle(grids);
    update_wirelength();
    ori_wire_length = wire_length;
}

void Net::update_wirelength() {
    unordered_set<Point, MyHashFunction> visited;
    wire_length = 0;
    for(auto& node : branch_nodes) {
        for(auto& neighbor : node.second.neighbors) {
            if(visited.find(neighbor.first) != visited.end())
                continue;
            neighbor.second.update_wire_length();
            wire_length+=neighbor.second.wire_length;
        }
        visited.insert(node.first);
    }
}

void Net::update_wirelength_fast() {
    unordered_set<Point, MyHashFunction> visited;
    wire_length = 0;
    for(auto& node : branch_nodes) {
        for(auto& neighbor : node.second.neighbors) {
            if(visited.find(neighbor.first) != visited.end())
                continue;
            wire_length+=neighbor.second.wire_length;
        }
        visited.insert(node.first);
    }
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
        std::vector<TwoPinNet>& routingTree,
        vector<vector<vector<Gcell>>>& grids
        ) {
    bool is_pin = (pin_map.find(start_p) != pin_map.end());
    while(true) {  // for 6 dir     
        Point dir = return_next_dir(degreeMap, start_p);
        if(dir == Point(-1,-1,-1))
            break;
        // -dir because delete outgoing edge
        add_net_demand_into_graph(start_p.x, start_p.y, start_p.z, grids);
        decrese_degree_endpoint(degreeMap, start_p, Point(-dir.x, -dir.y, -dir.z));
        Point now_p = start_p;
        TwoPinNet segment;
        segment.n1.p = start_p;
        if(segment.n1.mergedLocalPins.empty())
            segment.n1.mergedLocalPins = localNets[{now_p.x, now_p.y, now_p.z}];
        
        if(localNets[{now_p.x, now_p.y, now_p.z}].size() > 1) segment.n1.type = 2;
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
                add_net_demand_into_graph(now_p.x, now_p.y, now_p.z, grids);
                decrese_degree_endpoint(degreeMap, now_p, dir);
                path.second = now_p;
                segment.paths.push_back(path);
                segment.n2.p = now_p;
                if(segment.n2.mergedLocalPins.empty())
                    segment.n2.mergedLocalPins = localNets[{now_p.x, now_p.y, now_p.z}];
                if(localNets[{now_p.x, now_p.y, now_p.z}].size() > 1)
                    segment.n2.type = 2;
                else
                    segment.n2.type = (is_pin) ? 1 : 0;
                if(!(segment.n1.p == segment.n2.p))
                    routingTree.push_back(segment);
                else {
                    del_twoPinNet_from_graph(segment, grids);
                    if(!localNets[{segment.n1.p.x, segment.n1.p.y, segment.n1.p.z}].empty()) {
                        branch_nodes[segment.n1.p].node.mergedLocalPins = localNets[{segment.n1.p.x, segment.n1.p.y, segment.n1.p.z}];
                        if(branch_nodes[segment.n1.p].node.mergedLocalPins.size() > 1) branch_nodes[segment.n1.p].node.type = 2;
                        else branch_nodes[segment.n1.p].node.type = 1;
                    }
                }
                traverse_passing_map(degreeMap, pin_map, steiner_map, now_p, localNets, routingTree, grids);
                break;
            }
            else {
                // keep propogate
                Point next_p = now_p + dir;            
                if(check_map_legal(degreeMap, next_p) && check_map_dir(degreeMap, now_p, dir)) {
                    add_net_demand_into_graph(now_p.x, now_p.y, now_p.z, grids);
                    decrese_degree_middle_p(degreeMap, now_p, dir);
                }
                else {
                    // create segment path
                    add_net_demand_into_graph(now_p.x, now_p.y, now_p.z, grids);
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
                        add_net_demand_into_graph(now_p.x, now_p.y, now_p.z, grids);
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
        if(twopin.n1.type == 2 || twopin.n1.type == 1) {
            cout << "local net Node1: \n";
            for(auto& pin : twopin.n1.mergedLocalPins) cout << "Cell " << pin.first+1 << " Pin " << pin.second+1  << " ";
            cout << "\n";
        }
        if(twopin.n2.type == 2 || twopin.n2.type == 1) {
            cout << "local net Node2: \n";
            for(auto& pin : twopin.n2.mergedLocalPins) cout << "Cell " << pin.first+1 << " Pin " << pin.second+1 << " ";
            cout << "\n";
        }
        twopin.update_wire_length();
        cout << "wl: " << twopin.wire_length << endl;
    }  
}

void Net::print_branch_nodes() {
    for(auto& node : branch_nodes) {
        cout << node.first << " " << node.second.node.type << endl;
        for(auto& neighbor : node.second.neighbors) {
            cout << "\t" << neighbor.first << endl;
            for(auto& path : neighbor.second.paths)
                cout << "\t\t" << path.first << "->" << path.second << endl;
        }
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
        TwoPinNet r_twopin;
        r_twopin.n1 = twopin.n2;
        r_twopin.n2 = twopin.n1;
        r_twopin.paths.insert(r_twopin.paths.begin(), twopin.paths.rbegin(), twopin.paths.rend());
        for(auto& path : r_twopin.paths)
            swap(path.first, path.second);
        p2_ptr->second.neighbors.push_back(make_pair(p1,r_twopin));
    }  
}

void Net::remove_dangling_wire(vector<vector<vector<Gcell>>>& grids) {
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
        if(treenode.node.type == 0 && treenode.neighbors.size() == 2)
            clear_steiner_point(tree_p, grids);

        // find dangling endpoint
        if( treenode.node.type == -1 || (treenode.node.type == 0 && treenode.neighbors.size() == 1) ) {
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
            del_twoPinNet_from_graph(treenode.neighbors[0].second, grids);
            branch_nodes.erase(tree_p);
        }
    }
}

void Net::remove_branch_cycle(vector<vector<vector<Gcell>>>& grids) {
    std::priority_queue<TwoPinNet, vector<TwoPinNet>, greater<TwoPinNet>> frontier_edges;
    std::queue<TwoPinNet> skip_edges;
    std::unordered_set <Point,MyHashFunction> visited_nodes;
    unordered_set<Point, MyHashFunction> candidateSteiner;
    push_edge_in_queue(frontier_edges);
    // construct MST
    while(!frontier_edges.empty()) {
        auto edge = frontier_edges.top();
        frontier_edges.pop();
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
        if(e2_iter==visited_nodes.end() && e1_iter!=visited_nodes.end()) {
            visited_nodes.insert(edge.n2.p);  
            while(!skip_edges.empty()) {
                frontier_edges.push(skip_edges.front());
                skip_edges.pop();
            }
            continue;
        } 
        if(e1_iter!=visited_nodes.end() && e2_iter!=visited_nodes.end()) {
            //delete this edge in branch_nodes
            int wire_len = edge.wire_length;
            auto& e1_treeNode = branch_nodes[edge.n1.p];
            auto& e2_treeNode = branch_nodes[edge.n2.p];
            for(int n=0; n<e1_treeNode.neighbors.size(); n++) {
                e1_treeNode.neighbors[n].second.update_wire_length();
                if(e1_treeNode.neighbors[n].first == edge.n2.p && wire_len == e1_treeNode.neighbors[n].second.wire_length) {
                    del_twoPinNet_from_graph((e1_treeNode.neighbors.begin()+n)->second, grids); 
                    e1_treeNode.neighbors.erase(e1_treeNode.neighbors.begin()+n);
                    candidateSteiner.insert(edge.n1.p);
                    break;
                }
            }
            for(int n=0; n<e2_treeNode.neighbors.size(); n++) {
                e2_treeNode.neighbors[n].second.update_wire_length();
                if(e2_treeNode.neighbors[n].first == edge.n1.p && wire_len == e2_treeNode.neighbors[n].second.wire_length) {
                    e2_treeNode.neighbors.erase(e2_treeNode.neighbors.begin()+n);
                    candidateSteiner.insert(edge.n2.p);
                    break;
                }
            }
            continue;
        } else {
            skip_edges.push(edge);
        }
    }
    for(auto& p : candidateSteiner) {
        auto pos = branch_nodes.find(p);
        if(pos == branch_nodes.end()) continue;
        TreeNode& treeNode = pos->second;
        if(treeNode.node.type == 0 && treeNode.neighbors.size() == 2) clear_steiner_point(p, grids);
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
            frontier_edges.push(neighbor.second);
            traverse_points.insert(neighbor.first);
        }
        visited_points.insert(bfs_p);
    }
}

void Net::clear_steiner_point(Point p, vector<vector<vector<Gcell>>>& grids) {
    TreeNode& treeNode = branch_nodes[p];
    const Point& neighbor1 = treeNode.neighbors[0].first;
    const Point& neighbor2 = treeNode.neighbors[1].first;
    TwoPinNet& twoPinNet1 = treeNode.neighbors[0].second;
    TwoPinNet& twoPinNet2 = treeNode.neighbors[1].second;
    TreeNode& neighbor1TreeNode = branch_nodes[neighbor1];
    TreeNode& neighbor2TreeNode = branch_nodes[neighbor2];
    if(neighbor1 == neighbor2) {
        del_twoPinNet_from_graph(twoPinNet2, grids);
        for(int i = neighbor1TreeNode.neighbors.size() - 1; i >= 0; --i) {
            if(neighbor1TreeNode.neighbors[i].first == p) {
                neighbor1TreeNode.neighbors.erase(neighbor1TreeNode.neighbors.begin()+i);
                break;
            }
        }
        treeNode.neighbors.pop_back();
        treeNode.node.type = -1;
        return;
    }
    int pos1 = -1, pos2 = -1, pos3 = -1, pos4 = -1;
    for(int i = 0; i < neighbor1TreeNode.neighbors.size(); ++i) {
        if(neighbor1TreeNode.neighbors[i].first == p) pos1 = i;
        if(neighbor1TreeNode.neighbors[i].first == neighbor2) pos2 = i;
        if(pos1 != -1 && pos2 != -1) break;
    }
    for(int i = 0; i < neighbor2TreeNode.neighbors.size(); ++i) {
        if(neighbor2TreeNode.neighbors[i].first == p) pos3 = i;
        if(neighbor2TreeNode.neighbors[i].first == neighbor1) pos4 = i;
        if(pos3 != -1 && pos4 != -1) break;
    }
    //Multi Edge occurs
    if(pos2 != -1 && pos4 != -1 ) {
        int wirelength1, wirelength2;
        neighbor1TreeNode.neighbors[pos2].second.update_wire_length();
        wirelength1 = neighbor1TreeNode.neighbors[pos2].second.wire_length;
        twoPinNet1.update_wire_length();
        twoPinNet2.update_wire_length();
        wirelength2 = twoPinNet1.wire_length + twoPinNet2.wire_length;
        if(wirelength1 <= wirelength2) {
            del_twoPinNet_from_graph(twoPinNet1, grids);
            del_twoPinNet_from_graph(twoPinNet2, grids);
        }
        else {
            del_twoPinNet_from_graph(neighbor1TreeNode.neighbors[pos2].second, grids);
            //grids[p.x][p.y][p.z].passingNets[netId]--;
            neighbor1TreeNode.neighbors[pos2].second.paths = neighbor1TreeNode.neighbors[pos1].second.paths;
            neighbor2TreeNode.neighbors[pos4].second.paths = neighbor2TreeNode.neighbors[pos3].second.paths;
            neighbor1TreeNode.neighbors[pos2].second.paths.insert(neighbor1TreeNode.neighbors[pos2].second.paths.end(), twoPinNet2.paths.begin(), twoPinNet2.paths.end());
            neighbor2TreeNode.neighbors[pos4].second.paths.insert(neighbor2TreeNode.neighbors[pos4].second.paths.end(), twoPinNet1.paths.begin(), twoPinNet1.paths.end());
        }
    }
    else{
        //grids[p.x][p.y][p.z].passingNets[netId]--;
        neighbor1TreeNode.neighbors.push_back({neighbor2, TwoPinNet()});
        neighbor2TreeNode.neighbors.push_back({neighbor1, TwoPinNet()});
        int index1 = neighbor1TreeNode.neighbors.size() - 1;
        int index2 = neighbor2TreeNode.neighbors.size() - 1;
        neighbor1TreeNode.neighbors[index1].second.n1 = neighbor1TreeNode.node;
        neighbor1TreeNode.neighbors[index1].second.n2 = neighbor2TreeNode.node;
        neighbor2TreeNode.neighbors[index2].second.n1 = neighbor2TreeNode.node;
        neighbor2TreeNode.neighbors[index2].second.n2 = neighbor1TreeNode.node;
        neighbor1TreeNode.neighbors[index1].second.paths = neighbor1TreeNode.neighbors[pos1].second.paths;
        neighbor2TreeNode.neighbors[index2].second.paths = neighbor2TreeNode.neighbors[pos3].second.paths;
        neighbor1TreeNode.neighbors[index1].second.paths.insert(neighbor1TreeNode.neighbors[index1].second.paths.end(),
                                                                twoPinNet2.paths.begin(), twoPinNet2.paths.end());
        neighbor2TreeNode.neighbors[index2].second.paths.insert(neighbor2TreeNode.neighbors[index2].second.paths.end(),
                                                                twoPinNet1.paths.begin(), twoPinNet1.paths.end());
    }
    neighbor1TreeNode.neighbors.erase(neighbor1TreeNode.neighbors.begin()+pos1);
    neighbor2TreeNode.neighbors.erase(neighbor2TreeNode.neighbors.begin()+pos3);
    branch_nodes.erase(p);
}

bool Net::insert_steiner_point(Point p, vector<vector<vector<Gcell>>>& grids) {
    unordered_set<Point, MyHashFunction> used_points;
    for(auto& branch_node : this->branch_nodes) {
        used_points.insert(branch_node.first);
        for(auto& neighbor : branch_node.second.neighbors) {
            // traverse before
            if(used_points.find(neighbor.first) != used_points.end())
                continue;         
            vector<pair<Point,Point>>& ori_paths = neighbor.second.paths;
            for(int path_idx=0; path_idx<ori_paths.size(); path_idx++) {
                auto& path = ori_paths[path_idx];
                Point dir = norm(path.first - path.second);
                Point p_dir = (p == path.first) ? norm(p-path.second) : norm(p-path.first);
                Point tmp_first = path.first, tmp_second = path.second;
                if(dir != p_dir)
                    continue;
                if( (dir.x != 0 && in_range(p.x, path.first.x, path.second.x)) ||
                    (dir.y != 0 && in_range(p.y, path.first.y, path.second.y)) ||
                    (dir.z != 0 && in_range(p.z, path.first.z, path.second.z)) ||
                    (p == path.first) || (p == path.second) ) {                   
                    // divide path       
                    if(p != path.first && p != path.second) {               
                        path.second = p;
                        ori_paths.insert(ori_paths.begin()+path_idx+1, {p, tmp_second});
                        // add passing net at steiner point grid
                        grids[p.x][p.y][p.z].passingNets[netId]++;
                    }
                    vector<pair<Point,Point>> new_paths(ori_paths.begin()+path_idx+1, ori_paths.end());
                    ori_paths.erase(ori_paths.begin()+path_idx+1, ori_paths.end());
                    // insert new steiner point
                    Node new_steiner_p(p,0);
                    Point ori_left_p = branch_node.first, ori_right_p = neighbor.first;
                    branch_nodes[p].node = new_steiner_p;
                    branch_nodes[p].neighbors.resize(2);
                    branch_nodes[p].neighbors[0].first = ori_left_p;
                    branch_nodes[p].neighbors[1].first = ori_right_p;
                    branch_nodes[p].neighbors[0].second = neighbor.second;
                    Node tmp_right_node = neighbor.second.n2;
                    // find left neighbor to update
                    for(int n=0; n<branch_nodes[ori_left_p].neighbors.size(); n++) {
                        if(branch_nodes[ori_left_p].neighbors[n].first == ori_right_p) {
                            branch_nodes[ori_left_p].neighbors[n].first = p;
                            branch_nodes[ori_left_p].neighbors[n].second.n2.p = p;
                            branch_nodes[p].neighbors[0].second = two_pin_reverse(branch_nodes[ori_left_p].neighbors[n].second);
                        }
                    }
                    // update right two_pin
                    TwoPinNet& new_two_pin = branch_nodes[p].neighbors[1].second;
                    new_two_pin.n1 = new_steiner_p;
                    new_two_pin.n2 = tmp_right_node;
                    new_two_pin.paths = new_paths;
                    neighbor.second.n2 = new_steiner_p;
                    // find right neighbor to update
                    for(int n=0; n<branch_nodes[ori_right_p].neighbors.size(); n++) {
                        if(branch_nodes[ori_right_p].neighbors[n].first == ori_left_p) {
                            branch_nodes[ori_right_p].neighbors[n].first = p;
                            // reverse two pin
                            TwoPinNet r_twopin;
                            r_twopin.n1 = new_two_pin.n2;
                            r_twopin.n2 = new_two_pin.n1;
                            r_twopin.paths.insert(r_twopin.paths.begin(), new_two_pin.paths.rbegin(), new_two_pin.paths.rend());
                            for(auto& path : r_twopin.paths)
                                swap(path.first, path.second);
                            branch_nodes[ori_right_p].neighbors[n].second = r_twopin;
                        }
                    }
                    return 1;
                }
            }
        }
    }
    return 0; 
}

void Net::set_point_component(unordered_map<Point, int, MyHashFunction>& component_map) {
    int component_count = 0;
    for(auto& branch_node : branch_nodes) {
        if(component_map.find(branch_node.first) != component_map.end()) continue;
        component_map[branch_node.first] = ++component_count;
        queue<Point> nodes;
        nodes.push(branch_node.first);
        while(!nodes.empty()) {
            Point current_p = nodes.front();
            nodes.pop();
            for(auto& neighbor : branch_nodes[current_p].neighbors) {
                if(component_map.find(neighbor.first) == component_map.end()) {
                    component_map[neighbor.first] = component_count;
                    nodes.push(neighbor.first);
                }        
            }    
        }       
    }
}

RoutingGraph::RoutingGraph(): usedCellMove(0) {
    segmentTree = new SegmentTree(*this);
}

RoutingGraph::RoutingGraph(random_device& rd): usedCellMove(0) {
    gen = default_random_engine(rd());
    dis = uniform_real_distribution<>(0.0,1.0);
    segmentTree = new SegmentTree(*this);
}
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
    //Update merged local pin
    Cell& cell = cellInstances[cellIndex];
    if(!nets.empty()) {
        for(int i = 0; i < cell.pins.size(); ++i) {
            Pin& pin = cell.pins[i];
            int netIndex = pin.connectedNet;
            if(netIndex != -1) {
                Net& net = nets[netIndex];
                if(net.branch_nodes.find(Point(x, y, pin.layer)) == net.branch_nodes.end()) {
                    net.insert_steiner_point(Point(x, y, pin.layer), grids);
                }
                TreeNode& treeNode = net.branch_nodes[Point(x, y, pin.layer)];
                Node& node = treeNode.node;
                net.add_net_demand_into_graph(x, y, pin.layer, grids);
                node.mergedLocalPins.push_back({cellIndex, i});
                if(node.type == 0) node.type = 1;
                else if(node.type == 1 && node.mergedLocalPins.size() > 1) node.type = 2;
                if(pin.pseudo) {
                    //cout << "pseudo pin " << i << endl;
                    int actualLayer = pin.actualPinLayer;
                    for(int j = actualLayer; j <= pin.layer; ++j)
                        net.add_net_demand_into_graph(x, y, j, grids);
                }
            }
        }
    }
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
        if(netIndex != -1) {
            Net& net = nets[netIndex];
            net.del_net_from_graph(x, y, pin.layer, grids);
            net.branch_nodes.clear();
        }
    }
}

void RoutingGraph::del_cell_neighbor(int cellIndex) {
    Cell& cell = cellInstances[cellIndex];
    int x = cell.x;
    int y = cell.y;
    placement[x][y].erase(cellIndex);
    del_cell_demand_from_graph(x, y, cellInstances[cellIndex].mcType);
    for(int i = 0; i < cell.pins.size(); ++i) {
        Pin& pin = cell.pins[i];
        int netIndex = pin.connectedNet;
        Point p = Point(x,y,pin.layer);
        if(netIndex != -1) {
            Net& net = nets[netIndex];
            if(pin.pseudo) {
                int actualLayer = pin.actualPinLayer;
                for(int j = actualLayer; j <= pin.layer; ++j)
                    net.del_seg_demand_from_graph(x, y, j, grids);
            }
            if(net.branch_nodes.empty() || net.branch_nodes.find(p) == net.branch_nodes.end()) continue;

            TreeNode& treeNode = net.branch_nodes[p];

            for(auto it = treeNode.node.mergedLocalPins.begin(); it != treeNode.node.mergedLocalPins.end();) {
                if(it->first == cellIndex && it->second == i) {
                    it = treeNode.node.mergedLocalPins.erase(it);
                    net.del_seg_demand_from_graph(x, y, pin.layer, grids);
                }
                else ++it;
            }

            if(treeNode.node.mergedLocalPins.size() == 1)
                treeNode.node.type = 1;

            if(treeNode.node.mergedLocalPins.size() >= 1)
                continue;

            for(auto it = treeNode.neighbors.begin(); it != treeNode.neighbors.end(); ++it) {
                const Point& neighbor = it->first;
                net.wire_length -= it->second.wire_length;
                net.del_twoPinNet_from_graph(it->second, grids);
                TreeNode& neighborTreeNode = net.branch_nodes[neighbor];
                for(int j = 0; j < neighborTreeNode.neighbors.size(); ++j) {
                    if(neighborTreeNode.neighbors[j].first == p) {
                        neighborTreeNode.neighbors.erase(neighborTreeNode.neighbors.begin() + j);
                        break;
                    }
                }
            }
            net.branch_nodes.erase(p);
        }
    }
}

void RoutingGraph::del_cell_last_k_neighbor(int cellIndex, unordered_map<int, int>& netK) {
    Cell& cell = cellInstances[cellIndex];
    int x = cell.x;
    int y = cell.y;
    placement[x][y].erase(cellIndex);
    del_cell_demand_from_graph(x, y, cellInstances[cellIndex].mcType);
    for(int i = 0; i < cell.pins.size(); ++i) {
        Pin& pin = cell.pins[i];
        int netIndex = pin.connectedNet;
        Point p = Point(x,y,pin.layer);
        if(netIndex != -1) {
            Net& net = nets[netIndex];
            if(pin.pseudo) {
                int actualLayer = pin.actualPinLayer;
                for(int j = actualLayer; j <= pin.layer; ++j)
                    net.del_seg_demand_from_graph(x, y, j, grids);
            }
            if(net.branch_nodes.empty() || net.branch_nodes.find(p) == net.branch_nodes.end()) continue;
            TreeNode& treeNode = net.branch_nodes[p];
            for(auto it = treeNode.node.mergedLocalPins.begin(); it != treeNode.node.mergedLocalPins.end();) {
                if(it->first == cellIndex) {
                    it = treeNode.node.mergedLocalPins.erase(it);
                    net.del_seg_demand_from_graph(x, y, pin.layer, grids);
                }
                else ++it;
            }

            if(treeNode.node.mergedLocalPins.size() == 1)
                treeNode.node.type = 1;
            else if(treeNode.node.mergedLocalPins.size() == 0) treeNode.node.type = 0;

            if(treeNode.neighbors.empty() && treeNode.node.type == 0) {
                net.branch_nodes.erase(p);
                continue;
            }
            if(treeNode.neighbors.size() == 2 && treeNode.node.type == 0) {
                net.clear_steiner_point(p, grids);
                continue;
            }

            unordered_map<int, int>::iterator pos = netK.find(netIndex);
            if(pos == netK.end()) continue;
            int numberToDelete = pos->second;
            for(int j = 0; j < numberToDelete; j++) {
                pair<Point, TwoPinNet>& lastElement = treeNode.neighbors.back();
                Point& neighbor = lastElement.first;
                net.del_twoPinNet_from_graph(lastElement.second, grids);
                TreeNode& neighborTreeNode = net.branch_nodes[neighbor];
                for(int k = neighborTreeNode.neighbors.size() - 1; k >= 0; --k) {
                    if(neighborTreeNode.neighbors[k].first == p) {
                        neighborTreeNode.neighbors.erase(neighborTreeNode.neighbors.begin()+k);
                        break;
                    }
                }
                if(neighborTreeNode.node.type == 0 && neighborTreeNode.neighbors.size() == 2) {
                    net.clear_steiner_point(neighbor, grids);                   
                }
                treeNode.neighbors.pop_back();             
            }
            netK.erase(pos);
            if(treeNode.node.mergedLocalPins.size() >= 1)
                continue;

            if(treeNode.neighbors.empty())
                net.branch_nodes.erase(p);
            else
                treeNode.node.type = 0;           
        }
    }
}

void RoutingGraph::del_cell_point_net(int cellIndex, vector<pair<Point, int>>& point_nets) {
    Cell& cell = cellInstances[cellIndex];
    int x = cell.x;
    int y = cell.y;
    placement[x][y].erase(cellIndex);
    del_cell_demand_from_graph(x, y, cellInstances[cellIndex].mcType);
    // del designate net on point
    for(int i = point_nets.size()-1; i >= 0; i--) {
        pair<Point, int>& p_net = point_nets[i];
        Point p = p_net.first;
        int netIndex = p_net.second;
        Net& net = nets[netIndex];
        //cout << "\n#print seg\n";
        //net.print_branch_nodes();
        if(net.branch_nodes.find(p) != net.branch_nodes.end() && net.branch_nodes[p].neighbors.size()>0) {
            TreeNode& treeNode = net.branch_nodes[p];
            pair<Point, TwoPinNet>& lastElement = treeNode.neighbors.back();
            Point& neighbor = lastElement.first;
            net.del_twoPinNet_from_graph(lastElement.second, grids);
            TreeNode& neighborTreeNode = net.branch_nodes[neighbor];
            // find neighbor branch node to p, then delete
            for(int k = neighborTreeNode.neighbors.size() - 1; k >= 0; --k) {
                if(neighborTreeNode.neighbors[k].first == p) {
                    neighborTreeNode.neighbors.erase(neighborTreeNode.neighbors.begin()+k);
                    break;
                }
            }
            net.wire_length -= treeNode.neighbors.back().second.wire_length;
            treeNode.neighbors.pop_back();          
        }
    }
    // del cell
    for(int i = 0; i < cell.pins.size(); ++i) {
        Pin& pin = cell.pins[i];
        int netIndex = pin.connectedNet;
        Point p = Point(x,y,pin.layer);
        if(netIndex != -1) {
            Net& net = nets[netIndex];
            if(pin.pseudo) {
                int actualLayer = pin.actualPinLayer;
                for(int j = actualLayer; j <= pin.layer; ++j)
                    net.del_seg_demand_from_graph(x, y, j, grids);
            }
            if(net.branch_nodes.empty() || net.branch_nodes.find(p) == net.branch_nodes.end()) continue;
            TreeNode& treeNode = net.branch_nodes[p];
            for(auto it = treeNode.node.mergedLocalPins.begin(); it != treeNode.node.mergedLocalPins.end();) {
                if(it->first == cellIndex) {
                    it = treeNode.node.mergedLocalPins.erase(it);
                    net.del_seg_demand_from_graph(x, y, pin.layer, grids);
                }
                else ++it;
            }
            if(treeNode.node.mergedLocalPins.size() == 1)
                treeNode.node.type = 1;
            else if(treeNode.node.mergedLocalPins.size() == 0) treeNode.node.type = 0;
            if(treeNode.neighbors.empty() && treeNode.node.type == 0) {
                net.branch_nodes.erase(p);
                continue;
            }

            if(treeNode.node.mergedLocalPins.size() >= 1)
                continue;

            if(treeNode.neighbors.empty())
                net.branch_nodes.erase(p);
            else
                treeNode.node.type = 0;           
        }
    }
}

void Net::add_net_demand_into_graph(int x, int y, int z, vector<vector<vector<Gcell>>>& grids) {
    if(grids[x][y][z].passingNets[netId]++ == 0) //every pin and segment contribute 1 to passingNet
        grids[x][y][z].demand++;
}

void Net::add_twopin_demand_into_graph(TwoPinNet& twoPinNet, vector<vector<vector<Gcell>>>& grids) {
    for(auto it = twoPinNet.paths.begin(); it != twoPinNet.paths.end(); ++it) {
        Point& start = it->first;
        Point& end = it->second;
        if(start.x != end.x) {
            if(start.x > end.x) swap(start.x, end.x);
            for(int i = start.x; i <= end.x; ++i)
                add_net_demand_into_graph(i, start.y, start.z, grids);
        }
        else if(start.y != end.y) {
            if(start.y > end.y) swap(start.y, end.y);
            for(int i = start.y; i <= end.y; ++i)
                add_net_demand_into_graph(start.x, i, start.z, grids);
        }
        else if(start.z != end.z) {
            if(start.z > end.z) swap(start.z, end.z);
            for(int i = start.z; i <= end.z; ++i)
                add_net_demand_into_graph(start.x, start.y, i, grids);
        }
    }
}

void Net::del_net_from_graph(int x, int y, int z, vector<vector<vector<Gcell>>>& grids) {
    if(--grids[x][y][z].passingNets[netId] == 0) //every pin automatically contributes 1 to passingNets
        grids[x][y][z].demand--;
    else {
        auto it = branch_nodes.begin();
        while(it != branch_nodes.end()) {
            const Point& point = it->first;
            TreeNode& treeNode = it->second;
            for(auto _it = treeNode.neighbors.begin(); _it != treeNode.neighbors.end(); _it++) {
                Point& neighbor = _it->first;
                TwoPinNet& twoPinNet = _it->second;
                TreeNode& neighborTreeNode = branch_nodes[neighbor];
                del_twoPinNet_from_graph(twoPinNet, grids);
                //delete neighbor reverse edge
                for(auto _iit = neighborTreeNode.neighbors.begin(); _iit != neighborTreeNode.neighbors.end(); ++_iit) {
                    if(_iit->first == point) {
                        neighborTreeNode.neighbors.erase(_iit);
                        break;
                    }
                }
            }
            it = branch_nodes.erase(it);
        }
    }
}

void Net::del_seg_demand(std::pair<Point,Point> segment, vector<vector<vector<Gcell>>>& grids) {
    int startRow = segment.first.x, startColumn = segment.first.y, startLayer = segment.first.z;
    int endRow = segment.second.x, endColumn = segment.second.y, endLayer = segment.second.z;
    if(startRow > endRow) 
        swap(startRow, endRow);
    if(startColumn > endColumn)
        swap(startColumn, endColumn);
    if(startLayer > endLayer)
        swap(startLayer, endLayer);
    // handle vertical segment
    if(startRow != endRow) {
        for(int j = startRow; j <= endRow; j++) {
            del_seg_demand_from_graph(j, startColumn, startLayer, grids);
        }
    }
    // handle horizontal segment
    else if(startColumn != endColumn) {
        for(int j = startColumn; j <= endColumn; j++) {
            del_seg_demand_from_graph(startRow, j, startLayer, grids);
        }
    }
    // handle via segment
    else if(startLayer != endLayer) {
        for(int j = startLayer; j <= endLayer; j++) {
            del_seg_demand_from_graph(startRow, startColumn, j, grids);
        }
    }
}

void Net::del_seg_demand_from_graph(int x, int y, int z, vector<vector<vector<Gcell>>>& grids) {
    if(--grids[x][y][z].passingNets[netId] == 0)
        grids[x][y][z].demand--;
    if(grids[x][y][z].demand<0)
        cout << "\nFuck\n\n";
}

void Net::del_twoPinNet_from_graph(TwoPinNet& twoPinNet, vector<vector<vector<Gcell>>>& grids) {
    for(auto it = twoPinNet.paths.begin(); it != twoPinNet.paths.end(); ++it) {
        Point& start = it->first;
        Point& end = it->second;
        del_seg_demand({start, end}, grids);
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
        //cout << "\nNew net: " << a++ << "\n";
        net.convert_seg_to_2pin(degreeMap, cellInstances, masterCells, grids);
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
    return flutetree;
}

bool sortbysec(const pair<Point,int> &a, const pair<Point,int> &b) 
{ 
    return (a.second > b.second); 
} 

struct cmp_big{
    bool operator() (const double& lhs, const double& rhs) {return lhs > rhs;}
};

bool RoutingGraph::try_to_ripple_move_into_optimal(int cellIdx, int& net_wirelength, vector<double>& cellGains) {
    net_wirelength = 0;
    Cell& cell = cellInstances[cellIdx];
    vector<pair<Point, int>> cells_pos;
    bool opt_flag = find_optimal_pos_without_check_overflow(cell, cells_pos);
    if(!opt_flag) return false;
    int selectedCell;
    double maxGain = -10000;
    for(int i = 0; i < cells_pos.size(); ++i) {
        if(i == moveArea) break;
        unordered_set<int>& placedCells = placement[cells_pos[i].first.x][cells_pos[i].first.y];
        for(const int& idx : placedCells) {
            if(cellGains[idx] > maxGain) {
                maxGain = cellGains[idx];
                selectedCell = idx;
            }
        }
    }
    if(maxGain == -10000) return 0;
    Cell& rippleCell = cellInstances[selectedCell];
    int rippleOriX = rippleCell.x, rippleOriY = rippleCell.y;
    int ripple_wl = 0, wl = 0;
    vector<pair<Point, int>> target_cells_pos;
    opt_flag = find_force_pos(rippleCell, target_cells_pos);
    if(opt_flag) return false;
    bool ripple_success = 0, success = 0;
    for(int i = 0; i < target_cells_pos.size(); ++i) {
        if(i == moveArea) break;
        Point to_p(target_cells_pos[i].first.x, target_cells_pos[i].first.y, 0);
        ripple_success = move_cell_reroute_or_reverse(to_p, selectedCell, ripple_wl);
        if(ripple_success) break;
    }
    cout << maxGain << " " << ripple_success << endl;
    if(ripple_success) {
        for(int i = 0; i < cells_pos.size(); ++i) {
            Point to_p(cells_pos[i].first.x, cells_pos[i].first.y, 0);
            success = move_cell_reroute_or_reverse(to_p, cellIdx, wl);
            if(success) break;
        }
        if(success) {
            net_wirelength = wl + ripple_wl;
            return true;
        }
        else {
            Point oriP(rippleOriX, rippleOriY, 0);
            move_cell_reroute_or_reverse(oriP, selectedCell, ripple_wl);
        }
    }
    net_wirelength = 0;
    return false;
}

bool RoutingGraph::try_to_swap_into_force(int cellIdx, int& net_wirelength) {
    net_wirelength = 0;
    Cell& cell = cellInstances[cellIdx];
    vector<pair<Point, int>> cells_pos;
    bool opt_flag = find_force_pos_without_check_overflow(cell, cells_pos);
    if(!opt_flag) return false;
    //Point p(cells_pos[0].first.x, cells_pos[0].first.y, 0);
    //unordered_set<int>& targetCells = placement[p.x][p.y];
    unordered_set<int> targetCells;
    for(int i = 0; i < cells_pos.size(); i++) {
        if(i == swapArea) break;
        Point p(cells_pos[i].first.x, cells_pos[i].first.y, 0);
        unordered_set<int>& placedCells = placement[p.x][p.y];
        for(const int& idx : placedCells) targetCells.insert(idx);
    }
    map<int, unordered_set<int>, cmp_big> targetSwapPredictGain;
    for(auto& cellIdx : targetCells) {
        Cell& targetCell = cellInstances[cellIdx];
        vector<pair<Point, int>> target_cells_pos;
        bool targetOptFlag = find_optimal_pos_without_check_overflow(targetCell, target_cells_pos);
        if(!targetOptFlag) continue;
        //Point targetP(target_cells_pos[0].first.x, target_cells_pos[0].first.y, 0);
        for(int i = 0; i < target_cells_pos.size(); ++i) {
            if(i == swapArea) break;
            Point targetP(target_cells_pos[i].first.x, target_cells_pos[i].first.y, 0);
            int distance = (ADIFF(cell.x, targetP.x)) + (ADIFF(cell.y, targetP.y));
            targetSwapPredictGain[distance].insert(cellIdx);
        }
        //int distance = (ADIFF(cell.x, targetP.x)) + (ADIFF(cell.y, targetP.y));
        //targetSwapPredictGain[distance].insert(cellIdx);
    }
    if(targetSwapPredictGain.empty()) return false;
    unordered_set<int>& candidates = targetSwapPredictGain.begin()->second;
    int selectedCell, maxPin = 0;
    for(auto& idx : candidates) {
        Cell& candidate = cellInstances[idx];
        if(candidate.pins.size() > maxPin) {
            selectedCell = idx;
            maxPin = candidate.pins.size();
        }
    }
    return swap_two_cells(cellIdx, selectedCell, net_wirelength);
}

bool RoutingGraph::try_to_swap_into_optimal(int cellIdx, int& net_wirelength) {
    net_wirelength = 0;
    Cell& cell = cellInstances[cellIdx];
    vector<pair<Point, int>> cells_pos;
    bool opt_flag = find_optimal_pos_without_check_overflow(cell, cells_pos);
    if(!opt_flag) return false;
    //Point p(cells_pos[0].first.x, cells_pos[0].first.y, 0);
    //unordered_set<int>& targetCells = placement[p.x][p.y];
    unordered_set<int> targetCells;
    for(int i = 0; i < cells_pos.size(); i++) {
        if(i == swapArea) break;
        Point p(cells_pos[i].first.x, cells_pos[i].first.y, 0);
        unordered_set<int>& placedCells = placement[p.x][p.y];
        for(const int& idx : placedCells) targetCells.insert(idx);
    }
    map<int, unordered_set<int>, cmp_big> targetSwapPredictGain;
    for(auto& cellIdx : targetCells) {
        Cell& targetCell = cellInstances[cellIdx];
        vector<pair<Point, int>> target_cells_pos;
        bool targetOptFlag = find_optimal_pos_without_check_overflow(targetCell, target_cells_pos);
        if(!targetOptFlag) continue;
        //Point targetP(target_cells_pos[0].first.x, target_cells_pos[0].first.y, 0);
        for(int i = 0; i < target_cells_pos.size(); ++i) {
            if(i == swapArea) break;
            Point targetP(target_cells_pos[i].first.x, target_cells_pos[i].first.y, 0);
            int distance = (ADIFF(cell.x, targetP.x)) + (ADIFF(cell.y, targetP.y));
            targetSwapPredictGain[distance].insert(cellIdx);
        }
        //int distance = (ADIFF(cell.x, targetP.x)) + (ADIFF(cell.y, targetP.y));
        //targetSwapPredictGain[distance].insert(cellIdx);
    }
    if(targetSwapPredictGain.empty()) return false;
    unordered_set<int>& candidates = targetSwapPredictGain.begin()->second;
    int selectedCell, maxPin = 0;
    for(auto& idx : candidates) {
        Cell& candidate = cellInstances[idx];
        if(candidate.pins.size() > maxPin) {
            selectedCell = idx;
            maxPin = candidate.pins.size();
        }
    }
    return swap_two_cells(cellIdx, selectedCell, net_wirelength);
}

void RoutingGraph::wirelength_driven_move(int& wl_improve) {
    map<double, set<int>, cmp_big> bucket;
    vector<double> cellGains(cellInstances.size(), 0.0);
    vector<int> cellNets(cellInstances.size(), 0);
    for(auto& net : nets) {
        int length = 0;
        unordered_map<Point, bool, MyHashFunction> visited;
        stack<Point> candidate;
        if(!net.branch_nodes.empty()) {
            const Point& start = net.branch_nodes.begin()->first;
            candidate.push(start);
            visited[start] = true;
            while(!candidate.empty()) {
                Point& cur = candidate.top();
                candidate.pop();
                TreeNode& treeNode = net.branch_nodes[cur];
                for(auto it = treeNode.neighbors.begin(); it != treeNode.neighbors.end(); ++it) {
                    const Point& neighbor = it->first;
                    TwoPinNet& twoPinNet = it->second;
                    twoPinNet.update_wire_length();
                    length += twoPinNet.wire_length;
                    if(visited[neighbor] != true) {
                        visited[neighbor] = true;
                        candidate.push(neighbor);
                    }
                }
            }
        }
        for(auto it = net.pins.begin(); it != net.pins.end(); ++it) {
            int cellIdx = it->first;
            cellNets[cellIdx] += length;
        }
    }
    for(int i = 0; i < cellInstances.size(); i++) {
        Cell& cell = cellInstances[i];
        vector<pair<Point, int>> cells_pos;
        //bool opt_flag = find_optimal_pos(cell, cells_pos);
        bool opt_flag;
        if(case3) {
            opt_flag = find_optimal_pos(cell, cells_pos);
        }
        else {
            if(mode == 0 || mode == 2) opt_flag = find_optimal_pos(cell, cells_pos);
            else opt_flag = find_force_pos(cell, cells_pos);
        }
        int gain = 0, accLength = 0;
        if(opt_flag) {
            Point to_p(cells_pos[0].first.x, cells_pos[0].first.y, 0);
            unordered_map<int, bool> netCalculated;
            for(auto& pin : cell.pins) {
                int l = column, r = -1, u = -1, d = row, t = -1, b = layer;
                if(pin.connectedNet != -1 && netCalculated[pin.connectedNet] == false) {
                    Net& net = nets[pin.connectedNet];
                    netCalculated[net.netId] = true;
                    for(auto& ele : net.pins) {
                        if(ele.first == i) {
                            l = min(l, to_p.y);
                            r = max(r, to_p.y);
                            u = max(u, to_p.x);
                            d = min(d, to_p.x);
                            t = max(t, pin.layer);
                            b = min(b, pin.layer);
                        }
                        else {
                            l = min(l, cellInstances[ele.first].y);
                            r = max(r, cellInstances[ele.first].y);
                            u = max(u, cellInstances[ele.first].x);
                            d = min(d, cellInstances[ele.first].x);
                            t = max(t, cellInstances[ele.first].pins[ele.second].layer);
                            b = min(b, cellInstances[ele.first].pins[ele.second].layer);
                        }
                    }
                }
                accLength += (r-l) + (u-d) + (t-b);
            }
            //cout << i << " " << cellNets.size() << " " << cells_pos.size() << "\n";
            cellGains[i] += (cellNets[i] - accLength); //bigger better
            //cellGains[i] -= cell.fail_count*10;
            cellGains[i] -= cell.fail_count*failPenalty;
            //cellGains[i] -= (double)placement[cell.x][cell.y].size();
            //cellGains[i] += cells_pos[0].second/layerNormalize;
            for(int j = 0; j < cells_pos.size(); ++j) {
                if(j == predictRegionSize) break;
                //case1
                if(case1) cellGains[j] += cells_pos[j].second/layerNormalize;
                else cellGains[i] += cells_pos[j].second/layerNormalize;
            }
        }
    }

    for(int i = 0; i < cellInstances.size(); ++i) {
        //cout << i << " : " << cellGains[i] << endl;
        bucket[cellGains[i]].insert(i);
    }
    int count = 0;
    wl_improve = 0;
    while(!bucket.empty()) {
        //if(this->movedCell.size() >= maxCellMove) return;
        set<int>& candidateList = bucket.begin()->second;
        if(candidateList.empty()) {
            bucket.erase(bucket.begin());
            continue;
        }
        int candidate = *(candidateList.begin());
        candidateList.erase(candidate);
        cellGains[candidate] = -10000;
        if(candidateList.empty())
            bucket.erase(bucket.begin());
        int wl=0;
        bool success = move_cell_into_optimal_region(candidate, wl);
        wl_improve = success ? wl_improve + wl : wl_improve;
        if(success)
            count++;
        else if(mode == 0 && count < maxCellMove/(topK*2)) {
            wl = 0;
            success = try_to_swap_into_optimal(candidate, wl);
            wl_improve = success ? wl_improve + wl : wl_improve;
            if(success) count++;
        }
        /*else if(mode == 1 && count < maxCellMove/20) {
            wl = 0;
            success = try_to_swap_into_force(candidate, wl);
            wl_improve = success ? wl_improve + wl : wl_improve;
            if(success) count++;
        }*/
        if(mode != 2)
            if(count > maxCellMove/topK) return;
    }
}

bool RoutingGraph::move_cell_into_optimal_region(int cell_idx, int& net_wirelength) {
    //cout << "\n#cell " << cell_idx << endl;
    Cell& cell = cellInstances[cell_idx];
    if(movedCell.find(cell_idx) == movedCell.end() && this->movedCell.size() >= maxCellMove) return 0;
    if(!cell.movable) return 0;
    vector<pair<Point,int>> cells_pos;
    bool opt_flag = 0;
    if(mode == 2) {
        Point to_p(cell.x, cell.y, 0);
        return move_cell_reroute_or_reverse(to_p, cell_idx, net_wirelength);
    }
    if(mode == 0)
        opt_flag = find_optimal_pos(cell, cells_pos);
    else if(mode == 1)
        opt_flag = find_force_pos(cell, cells_pos);
    if(cells_pos.size() == 0)
        return 0;
    if(!opt_flag)
        return 0;
    bool success = 0;
    for(int c=0; c<cells_pos.size(); c++) {
        if(c > moveArea) break;
        Point to_p(cells_pos[c].first.x,cells_pos[c].first.y,0);
        success = move_cell_reroute_or_reverse(to_p, cell_idx, net_wirelength);
        if(success) break;
    }
    cell.fail_count = (!success) ? cell.fail_count+1 : cell.fail_count;
    if(!success) {
        Point to_p(cell.x, cell.y, 0);
        move_cell_reroute_or_reverse(to_p, cell_idx, net_wirelength);
    }
    return success;
}

bool RoutingGraph::move_cell_reroute_or_reverse(Point to_p, int cell_idx, int& net_wirelength) {
    Cell& cell = cellInstances[cell_idx];
    int cell_ori_x = cell.x, cell_ori_y = cell.y;
    // key is net, value is tuple
    // source, sink, netId, source.node, twopin
    unordered_map<int, vector<tuple<Point,Point,int,Node,TwoPinNet>> > open_nets;
    // key is netId, value is branch_node
    unordered_map<int,unordered_map<Point, TreeNode, MyHashFunction>> branchs_copy;
    unordered_set<int> updated_branchs;
    unordered_map<int, map<tuple<int, int, int, int, int, int>, bool>> alreadyInOpenNets;
    unordered_set<int> net_id_pool;
    for(auto& pin : cell.pins) {     
        if(pin.connectedNet == -1)
            continue;
        auto& net = nets[pin.connectedNet];
        net_id_pool.insert(pin.connectedNet);
        // copy branch_node for reverse
        if(branchs_copy.find(pin.connectedNet) == branchs_copy.end())
            branchs_copy[pin.connectedNet] = net.branch_nodes;
        Point cell_p(cell.x, cell.y, pin.layer);
        unordered_set<Point, MyHashFunction> used_local_node;
        for(auto& neighbor : net.branch_nodes[cell_p].neighbors) {
            Point neighbor_p = neighbor.first;
            // local pin can not put in open net twice
            if(used_local_node.find(neighbor_p) != used_local_node.end()) {
                //cout << "skip\n";
                continue;
            }
            if(net.branch_nodes[neighbor_p].node.mergedLocalPins.size() >= 2)
                used_local_node.insert(neighbor_p);
            if(alreadyInOpenNets[net.netId][{to_p.x,to_p.y,pin.layer,neighbor_p.x, neighbor_p.y, neighbor_p.z}] == true) continue;
            neighbor.second.update_wire_length();
            net_wirelength += neighbor.second.wire_length;
            open_nets[net.netId].emplace_back(Point(to_p.x,to_p.y,pin.layer), neighbor_p, net.netId, net.branch_nodes[cell_p].node, neighbor.second);
            alreadyInOpenNets[net.netId][{to_p.x,to_p.y,pin.layer,neighbor_p.x, neighbor_p.y, neighbor_p.z}] = true;
        }        
        if(net.branch_nodes[cell_p].node.type == 2 && net.branch_nodes[cell_p].neighbors.empty()
        && alreadyInOpenNets[net.netId][{to_p.x,to_p.y,pin.layer,cell_ori_x, cell_ori_y, pin.layer}] == false) {
            open_nets[net.netId].emplace_back(Point(to_p.x,to_p.y,pin.layer), Point(cell_ori_x,cell_ori_y,pin.layer), net.netId, Node(), TwoPinNet());
            alreadyInOpenNets[net.netId][{to_p.x,to_p.y,pin.layer,cell_ori_x, cell_ori_y, pin.layer}] = true;
        }
    }
    int all_net_wl = 0;
    for(auto& netId : net_id_pool) {
        nets[netId].update_wirelength();
        all_net_wl += nets[netId].wire_length;
    }
    // delete cell neighbor two-pin
    del_cell_neighbor(cell_idx);
    add_cell(to_p.x,to_p.y,cell_idx);
    //cout << "#move " << cell_ori_x << " " << cell_ori_y << " to " << to_p.x << " " << to_p.y << endl;
    // test reroute
    vector<pair<Point, int>> point_nets;
    unordered_map<int, int> netK;
    bool routing_success = connect_all_nets(open_nets, net_wirelength, netK, point_nets);
    if(routing_success) {
        for(auto& netId : net_id_pool) {
            nets[netId].update_wirelength();
            all_net_wl -= nets[netId].wire_length;
        }
        if(openSA) {
            if(all_net_wl <= 0 && mode == 2) routing_success = 0;
            else if(all_net_wl <= 0) {
                double exponent = (double)(all_net_wl*1000-100)/temperature;
                exponent = exp(exponent);
                double prob = dis(gen);
                //cout << prob << " " << exponent << "\n";
                if(prob > exponent) {
                    routing_success = 0;
                }
                else net_wirelength = all_net_wl;
            }
            else net_wirelength = all_net_wl;
        }
        else {
            if(all_net_wl <= 0) routing_success = 0;
            else net_wirelength = all_net_wl;
        }
    }
    if(routing_success) {
        if(cell.x != cell.originalX || cell.y != cell_ori_y)
            movedCell.insert(cell_idx);
        return 1;
    }
    else {
        // reverse original net
        net_wirelength = 0;
        del_cell_point_net(cell_idx, point_nets);
        if(cell_ori_x == cell.originalX && cell_ori_y == cell.originalY)
            movedCell.erase(cell_idx);
        add_cell(cell_ori_x,cell_ori_y,cell_idx);
        for(auto& net_open_net : open_nets) {
            for(auto& open_net : net_open_net.second) {
                auto& net = nets[get<2>(open_net)];
                // resume net.branch_nodes
                if(updated_branchs.find(net.netId) == updated_branchs.end()) {
                    updated_branchs.insert(net.netId);
                    net.branch_nodes = branchs_copy[net.netId];
                }
                Point p1(cell_ori_x, cell_ori_y, get<0>(open_net).z);
                if(p1 == get<1>(open_net))
                    continue;
                // add two_pin demand into graph
                net.add_twopin_demand_into_graph(get<4>(open_net), grids);
                net.wire_length += get<4>(open_net).wire_length;
            }
        }
        //cout << "#after end\n";
        return 0;
    }
}

TwoPinNet RoutingGraph::convert_path_to_twopin(Point source, Point sink, unordered_map<Point,Point,MyHashFunction>& visited_p) {
    Point tmp_node = sink;       
    int dir = 0; // 1: vertical, 2: horizontal, 3:via

    TwoPinNet two_pin;
    two_pin.n1.p = source;
    two_pin.n2.p = sink;
    two_pin.paths.push_back({sink,visited_p[sink]});
    // build TwopinNet
    while(true) {
        if(tmp_node == source)
            break;
        tmp_node = visited_p[tmp_node];  
    }
    tmp_node = sink;   
    while(true) {
        Point back_p = visited_p[tmp_node];
        int new_dir = 0;
        if(back_p.x != tmp_node.x) new_dir = 1;
        if(back_p.y != tmp_node.y) new_dir = 2;
        if(back_p.z != tmp_node.z) new_dir = 3;
        if(new_dir != dir) {
            two_pin.paths.back().second = tmp_node;
            if(tmp_node!=source && tmp_node != sink)
                two_pin.paths.push_back({tmp_node,back_p});
            dir = new_dir;
        }
        if(tmp_node == source)
            break;
        tmp_node = back_p;
    }
    two_pin.update_wire_length();
    return two_pin_reverse(two_pin);
}

TwoPinNet RoutingGraph::convert_path_to_twopin_t2t(Point& source, Point sink, Point b_min, unordered_map<Point,Point,MyHashFunction>& visited_p,
    unordered_set<int>& source_comp_set, vector<vector<vector<int>>>& comp_grid_map) {
    Point tmp_node = sink;       
    int dir = 0; // 1: vertical, 2: horizontal, 3:via
    TwoPinNet two_pin;
    two_pin.n2.p = sink;
    two_pin.paths.push_back({sink,visited_p[sink]});
    while(true) {
        Point back_p = visited_p[tmp_node];
        //cout << "back_p: " << back_p << endl;
        int new_dir = 0;
        if(back_p.x != tmp_node.x) new_dir = 1;
        if(back_p.y != tmp_node.y) new_dir = 2;
        if(back_p.z != tmp_node.z) new_dir = 3;
        int tmp_comp = comp_grid_map[tmp_node.x-b_min.x][tmp_node.y-b_min.y][tmp_node.z-b_min.z];
        //cout << "tmp: " << tmp_node << " " << tmp_comp << endl;
        if(source_comp_set.find(tmp_comp) != source_comp_set.end()) {
            // find source
            //cout << "Find\n";
            source = tmp_node;
            two_pin.n1.p = tmp_node;
            two_pin.paths.back().second = source;
            break;
        }
        if(new_dir != dir) {
            two_pin.paths.back().second = tmp_node;
            if(tmp_node != sink)
                two_pin.paths.push_back({tmp_node,back_p});
            dir = new_dir;
        }
        tmp_node = back_p;
    }
    two_pin.update_wire_length();
    return two_pin_reverse(two_pin);
}

bool RoutingGraph::find_optimal_pos(Cell& cell, vector<pair<Point,int>>& cells_pos) {
    if(cell.movable == 0)
        return 0;
    // move cell to free space
    vector<int> x_series, y_series;
    for(auto& pin : cell.pins) {
        // find optimal region
        if(pin.connectedNet == -1)
            continue;
        Net& neighbor_net = nets[pin.connectedNet];     
        int min_x=INT_MAX, max_x=0, min_y=INT_MAX, max_y=0;    
        bool same_flag = 0;   
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
    if(x_series.size()==0 || y_series.size()==0) return 0;
    sort(x_series.begin(), x_series.end());
    sort(y_series.begin(), y_series.end());
    // count median boundary
    int opt_x_left = (x_series.size()%2) ? x_series[x_series.size()/2] : x_series[(x_series.size()-1)/2];
    int opt_x_right = (x_series.size()%2) ? x_series[x_series.size()/2] : x_series[(x_series.size()-1)/2+1];
    int opt_y_left = (y_series.size()%2) ? y_series[y_series.size()/2] : y_series[(y_series.size()-1)/2];
    int opt_y_right = (y_series.size()%2) ? y_series[y_series.size()/2] : y_series[(y_series.size()-1)/2+1];
    // cell already in optimal region
    if(cell.x >= opt_x_left && cell.x <= opt_x_right && cell.y >= opt_y_left && cell.y <= opt_y_right)
        return 0;
    for(int x=opt_x_left; x<=opt_x_right; x++) {
        for(int y=opt_y_left; y<=opt_y_right; y++) {
            int profit = check_cell_cost_in_graph(x, y, cell);
            if(profit > 0)
                cells_pos.emplace_back(Point(x,y,0), profit);
        }
    }
    if(cells_pos.empty()) return 0;
    sort(cells_pos.begin(), cells_pos.end(), sortbysec);
    return 1;
}

bool RoutingGraph::find_optimal_pos_without_check_overflow(Cell& cell, vector<pair<Point, int>>& cells_pos) {
    if(cell.movable == 0)
        return 0;
    // move cell to free space
    vector<int> x_series, y_series;
    for(auto& pin : cell.pins) {
        // find optimal region
        if(pin.connectedNet == -1)
            continue;
        Net& neighbor_net = nets[pin.connectedNet];
        int min_x=INT_MAX, max_x=0, min_y=INT_MAX, max_y=0;
        bool same_flag = 0;
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
    if(x_series.size()==0 || y_series.size()==0) return 0;
    sort(x_series.begin(), x_series.end());
    sort(y_series.begin(), y_series.end());
    // count median boundary
    int opt_x_left = (x_series.size()%2) ? x_series[x_series.size()/2] : x_series[(x_series.size()-1)/2];
    int opt_x_right = (x_series.size()%2) ? x_series[x_series.size()/2] : x_series[(x_series.size()-1)/2+1];
    int opt_y_left = (y_series.size()%2) ? y_series[y_series.size()/2] : y_series[(y_series.size()-1)/2];
    int opt_y_right = (y_series.size()%2) ? y_series[y_series.size()/2] : y_series[(y_series.size()-1)/2+1];
    // cell already in optimal region
    if(cell.x >= opt_x_left && cell.x <= opt_x_right && cell.y >= opt_y_left && cell.y <= opt_y_right)
        return 0;
    for(int x=opt_x_left; x<=opt_x_right; x++) {
        for(int y=opt_y_left; y<=opt_y_right; y++) {
            int profit = calculate_profit_without_check_overflow(x, y, cell);
            cells_pos.emplace_back(Point(x,y,0), profit);
        }
    }
    if(cells_pos.empty()) return 0;
    sort(cells_pos.begin(), cells_pos.end(), sortbysec);
    return 1;
}

bool RoutingGraph::find_force_pos(Cell& cell, vector<pair<Point,int>>& cells_pos) {
    if(cell.movable == 0)
        return 0;
    // move cell to free space
    vector<int> x_series, y_series;
    for(auto& pin : cell.pins) {
        // find optimal region
        if(pin.connectedNet == -1)
            continue;
        Net& net = nets[pin.connectedNet];
        Point p(cell.x, cell.y, pin.layer);
        if(net.branch_nodes.find(p) == net.branch_nodes.end())
            continue;
        for(auto& neighbor : net.branch_nodes[p].neighbors) {
            x_series.push_back(neighbor.first.x);
            y_series.push_back(neighbor.first.y);
        }
    }
    // no optimal region
    if(x_series.size()==0 || y_series.size()==0) return 0;
    sort(x_series.begin(), x_series.end());
    sort(y_series.begin(), y_series.end());
    // count median boundary
    int opt_x_left = (x_series.size()%2) ? x_series[x_series.size()/2] : x_series[(x_series.size()-1)/2];
    int opt_x_right = (x_series.size()%2) ? x_series[x_series.size()/2] : x_series[(x_series.size()-1)/2+1];
    int opt_y_left = (y_series.size()%2) ? y_series[y_series.size()/2] : y_series[(y_series.size()-1)/2];
    int opt_y_right = (y_series.size()%2) ? y_series[y_series.size()/2] : y_series[(y_series.size()-1)/2+1];
    // cell already in optimal region
    if(cell.x >= opt_x_left && cell.x <= opt_x_right && cell.y >= opt_y_left && cell.y <= opt_y_right)
        return 0;
    for(int x=opt_x_left; x<=opt_x_right; x++) {
        for(int y=opt_y_left; y<=opt_y_right; y++) {
            int profit = check_cell_cost_in_graph(x, y, cell);
            if(profit > 0)
                cells_pos.emplace_back(Point(x,y,0), profit);
        }
    }
    if(cells_pos.empty()) return 0;
    sort(cells_pos.begin(), cells_pos.end(), sortbysec);
    return 1;
}

bool RoutingGraph::find_force_pos_without_check_overflow(Cell& cell, vector<pair<Point, int>>& cells_pos) {
    if(cell.movable == 0)
        return 0;
    // move cell to free space
    vector<int> x_series, y_series;
    for(auto& pin : cell.pins) {
        // find optimal region
        if(pin.connectedNet == -1)
            continue;
        Net& net = nets[pin.connectedNet];
        Point p(cell.x, cell.y, pin.layer);
        if(net.branch_nodes.find(p) == net.branch_nodes.end())
            continue;
        for(auto& neighbor : net.branch_nodes[p].neighbors) {
            x_series.push_back(neighbor.first.x);
            y_series.push_back(neighbor.first.y);
        }
    }
    // no optimal region
    if(x_series.size()==0 || y_series.size()==0) return 0;
    sort(x_series.begin(), x_series.end());
    sort(y_series.begin(), y_series.end());
    // count median boundary
    int opt_x_left = (x_series.size()%2) ? x_series[x_series.size()/2] : x_series[(x_series.size()-1)/2];
    int opt_x_right = (x_series.size()%2) ? x_series[x_series.size()/2] : x_series[(x_series.size()-1)/2+1];
    int opt_y_left = (y_series.size()%2) ? y_series[y_series.size()/2] : y_series[(y_series.size()-1)/2];
    int opt_y_right = (y_series.size()%2) ? y_series[y_series.size()/2] : y_series[(y_series.size()-1)/2+1];
    // cell already in optimal region
    if(cell.x >= opt_x_left && cell.x <= opt_x_right && cell.y >= opt_y_left && cell.y <= opt_y_right)
        return 0;
    for(int x=opt_x_left; x<=opt_x_right; x++) {
        for(int y=opt_y_left; y<=opt_y_right; y++) {
            int profit = calculate_profit_without_check_overflow(x, y, cell);
            if(profit > 0)
                cells_pos.emplace_back(Point(x,y,0), profit);
        }
    }
    if(cells_pos.empty()) return 0;
    sort(cells_pos.begin(), cells_pos.end(), sortbysec);
    return 1;

}

int RoutingGraph::check_cell_cost_in_graph(int x, int y, Cell& cell) {
    int MCtype = cell.mcType;
    vector<int> layer_remain(layer), passing_net_profit(layer,0), pre_layer_remain(layer, 0), nxt_layer_remain(layer, 0);
    vector<unordered_set<int>> numberOfNetsOnLayer(layer);
    MasterCell& masterCell = masterCells[MCtype];
    for(auto& pin : cell.pins) {
        int netIndex = pin.connectedNet;
        if(netIndex == -1) continue;
        if(pin.pseudo) {
            for(int i = pin.actualPinLayer; i <= pin.layer; ++i) {
                auto pos = grids[x][y][i].passingNets.find(netIndex);
                if(pos == grids[x][y][i].passingNets.end() || pos->second == 0)
                    numberOfNetsOnLayer[i].insert(netIndex);
            }
        }
        else {
            auto pos = grids[x][y][pin.layer].passingNets.find(netIndex);
            if(pos == grids[x][y][pin.layer].passingNets.end() || pos->second == 0)
                numberOfNetsOnLayer[pin.layer].insert(netIndex);
        }
    }

    // initial remain
    for(int n=0; n<layer_remain.size(); n++) {
        layer_remain[n] = grids[x][y][n].capacity - grids[x][y][n].demand - numberOfNetsOnLayer[n].size();
        if(y > 0) pre_layer_remain[n] = grids[x][y-1][n].capacity - grids[x][y-1][n].demand;
        if(y < (column-1)) nxt_layer_remain[n] = grids[x][y+1][n].capacity - grids[x][y+1][n].demand;
    }
    // calculate locate at passing net benifit
    for(auto& pin : cell.pins) {
        int& p_layer = pin.layer;
        int& p_net_id = pin.connectedNet;
        if(grids[x][y][p_layer].passingNets.find(p_net_id) != grids[x][y][p_layer].passingNets.end())
            passing_net_profit[p_layer] += 50;
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
        }
        // check right cell demand
        if(y < (column - 1)) {
            nxtOriginalPairCnt = min(cellCount[x][y+1][MCtype2], cellCount[x][y][MCtype]);
            nxtAfterPairCnt = min(cellCount[x][y+1][MCtype2], cellCount[x][y][MCtype] + 1);
        }
        // Update current
        for(auto _it = it->second.begin(); _it != it->second.end(); ++_it) {
            pre_layer_remain[_it->first] -= _it->second * (preAfterPairCnt - preOriginalPairCnt);
            layer_remain[_it->first] -= _it->second * (preAfterPairCnt - preOriginalPairCnt + nxtAfterPairCnt - nxtOriginalPairCnt);
            nxt_layer_remain[_it->first] -= _it->second * (nxtAfterPairCnt - nxtOriginalPairCnt);
        }
    }
    // calculate profit
    int profit = 0;
    for(int n=0; n<layer_remain.size(); n++) {
        if(layer_remain[n] < 0 || pre_layer_remain[n] < 0 || nxt_layer_remain[n] < 0)
            return -1;
        profit += (layer_remain[n] + passing_net_profit[n]);
    }
    return profit;
}

int RoutingGraph::calculate_profit_without_check_overflow(int x, int y, Cell& cell) {
    int MCtype = cell.mcType;
    vector<int> layer_remain(layer), passing_net_profit(layer,0), pre_layer_remain(layer, 0), nxt_layer_remain(layer, 0);
    vector<unordered_set<int>> numberOfNetsOnLayer(layer);
    MasterCell& masterCell = masterCells[MCtype];
    for(auto& pin : cell.pins) {
        int netIndex = pin.connectedNet;
        if(netIndex == -1) continue;
        if(pin.pseudo) {
            for(int i = pin.actualPinLayer; i <= pin.layer; ++i) {
                auto pos = grids[x][y][i].passingNets.find(netIndex);
                if(pos == grids[x][y][i].passingNets.end() || pos->second == 0)
                    numberOfNetsOnLayer[i].insert(netIndex);
            }
        }
        else {
            auto pos = grids[x][y][pin.layer].passingNets.find(netIndex);
            if(pos == grids[x][y][pin.layer].passingNets.end() || pos->second == 0)
                numberOfNetsOnLayer[pin.layer].insert(netIndex);
        }
    }

    // initial remain
    for(int n=0; n<layer_remain.size(); n++) {
        layer_remain[n] = grids[x][y][n].capacity - grids[x][y][n].demand - numberOfNetsOnLayer[n].size();
        if(y > 0) pre_layer_remain[n] = grids[x][y-1][n].capacity - grids[x][y-1][n].demand;
        if(y < (column-1)) nxt_layer_remain[n] = grids[x][y+1][n].capacity - grids[x][y+1][n].demand;
    }
    // calculate locate at passing net benifit
    for(auto& pin : cell.pins) {
        int& p_layer = pin.layer;
        int& p_net_id = pin.connectedNet;
        if(grids[x][y][p_layer].passingNets.find(p_net_id) != grids[x][y][p_layer].passingNets.end())
            passing_net_profit[p_layer] += 50;
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
        }
        // check right cell demand
        if(y < (column - 1)) {
            nxtOriginalPairCnt = min(cellCount[x][y+1][MCtype2], cellCount[x][y][MCtype]);
            nxtAfterPairCnt = min(cellCount[x][y+1][MCtype2], cellCount[x][y][MCtype] + 1);
        }
        // Update current
        for(auto _it = it->second.begin(); _it != it->second.end(); ++_it) {
            pre_layer_remain[_it->first] -= _it->second * (preAfterPairCnt - preOriginalPairCnt);
            layer_remain[_it->first] -= _it->second * (preAfterPairCnt - preOriginalPairCnt + nxtAfterPairCnt - nxtOriginalPairCnt);
            nxt_layer_remain[_it->first] -= _it->second * (nxtAfterPairCnt - nxtOriginalPairCnt);
        }
    }
    // calculate profit
    int profit = 0;
    bool flag = false;
    for(int n=0; n<layer_remain.size(); n++) {
        if(layer_remain[n] < 0 || pre_layer_remain[n] < 0 || nxt_layer_remain[n] < 0)
            flag = true;
        profit += (layer_remain[n] + passing_net_profit[n]);
    }
    if(flag && profit > 0) return -1;
    return profit;
}

bool RoutingGraph::A_star_routing(Point source, Point sink, int NetId, unordered_map<Point,Point,MyHashFunction>& visited_p) {
    int min_x = min(source.x, sink.x);
    int max_x = max(source.x, sink.x);
    int min_y = min(source.y, sink.y);
    int max_y = max(source.y, sink.y);
    int min_l = min(source.z, sink.z);
    int max_l = max(source.z, sink.z);
    visited_p[source] = source;
    if(source == sink) return 0;
    max_l = (max_l<this->layer-1) ? max_l+1 : max_l;
    int min_r_l = nets[NetId].minRoutingLayer;
    //if(min_r_l>0) return 0;
    Point ori_source = source, ori_sink = sink;
    if(min_l < min_r_l) {
        min_l = min_r_l;
        if(source.z < min_r_l) {
            source.z = min_r_l;
            visited_p[source] = ori_source;
        }
        if(sink.z < min_r_l) {
            sink.z = min_r_l;
            visited_p[ori_sink] = sink;
        }
    }
    if(max_l < min_l) return 0;
    priority_queue<pair<Point,int>, vector<pair<Point,int>>, greater<pair<Point,int>>> p_q;
    Gcell& gcell = grids[source.x][source.y][source.z];
    int wire_length = 1;
    auto pos = gcell.passingNets.find(NetId);
    if(pos != gcell.passingNets.end() && pos->second != 0)
        wire_length = 0;
    int remain = gcell.capacity-gcell.demand-wire_length;
    if(remain < 0)
        return 0;
    p_q.emplace(source, gcell.demand+wire_length+distance(source,sink));
    bool find_flag = 0;
    while(!p_q.empty()) {
        auto frontier = p_q.top();
        p_q.pop();
        auto& f_point = frontier.first;
        if(f_point == sink) {
            find_flag = 1;
            break;
        }
        // up
        if(f_point.z%2==1 && f_point.x < max_x) {
            Point new_p(f_point.x+1,f_point.y,f_point.z);
            Gcell& new_gcell = grids[f_point.x+1][f_point.y][f_point.z];
            int wire_length = 1;
            auto new_pos = new_gcell.passingNets.find(NetId);
            if(new_pos != new_gcell.passingNets.end() && new_pos->second != 0)
                wire_length = 0;
            int remain = new_gcell.capacity-new_gcell.demand-wire_length;
            if(remain > 0 && visited_p.find(new_p) == visited_p.end()) {              
                p_q.emplace(new_p, new_gcell.demand+wire_length+distance(new_p,sink));
                visited_p[new_p] = f_point;
            }
        }
        // down
        if(f_point.z%2==1 && f_point.x > min_x) {
            Point new_p(f_point.x-1,f_point.y,f_point.z);
            Gcell& new_gcell = grids[f_point.x-1][f_point.y][f_point.z];
            int wire_length = 1;
            auto new_pos = new_gcell.passingNets.find(NetId);
            if(new_pos != new_gcell.passingNets.end() && new_pos->second != 0)
                wire_length = 0;
            int remain = new_gcell.capacity-new_gcell.demand-wire_length;
            if(remain > 0 && visited_p.find(new_p) == visited_p.end()) {
                p_q.emplace(new_p, new_gcell.demand+wire_length+distance(new_p,sink));
                visited_p[new_p] = f_point;
            }
        }
        // right
        if(f_point.z%2==0 && f_point.y < max_y) {
            Point new_p(f_point.x,f_point.y+1,f_point.z);
            Gcell& new_gcell = grids[f_point.x][f_point.y+1][f_point.z];
            int wire_length = 1;
            auto new_pos = new_gcell.passingNets.find(NetId);
            if(new_pos != new_gcell.passingNets.end() && new_pos->second != 0)
                wire_length = 0;
            int remain = new_gcell.capacity-new_gcell.demand-wire_length;
            if(remain > 0 && visited_p.find(new_p) == visited_p.end()) {
                p_q.emplace(new_p, new_gcell.demand+wire_length+distance(new_p,sink));
                visited_p[new_p] = f_point;
            }
        }
        // left
        if(f_point.z%2==0 && f_point.y > min_y) {
            Point new_p(f_point.x,f_point.y-1,f_point.z);
            Gcell& new_gcell = grids[f_point.x][f_point.y-1][f_point.z];
            int wire_length = 1;
            auto new_pos = new_gcell.passingNets.find(NetId);
            if(new_pos != new_gcell.passingNets.end() && new_pos->second != 0)
                wire_length = 0;
            int remain = new_gcell.capacity-new_gcell.demand-wire_length;
            if(remain > 0 && visited_p.find(new_p) == visited_p.end()) {
                p_q.emplace(new_p, new_gcell.demand+wire_length+distance(new_p,sink));
                visited_p[new_p] = f_point;
            }
        }
        // top
        if(f_point.z < max_l) {
            Point new_p(f_point.x,f_point.y,f_point.z+1);
            Gcell& new_gcell = grids[f_point.x][f_point.y][f_point.z+1];
            int wire_length = 1;
            auto new_pos = new_gcell.passingNets.find(NetId);
            if(new_pos != new_gcell.passingNets.end() && new_pos->second != 0)
                wire_length = 0;
            int remain = new_gcell.capacity-new_gcell.demand-wire_length;
            if(remain > 0 && visited_p.find(new_p) == visited_p.end()) {
                p_q.emplace(new_p, new_gcell.demand+wire_length+distance(new_p,sink));
                visited_p[new_p] = f_point;
            }
        }
        // down
        if(f_point.z > min_l) {
            Point new_p(f_point.x,f_point.y,f_point.z-1);
            Gcell& new_gcell = grids[f_point.x][f_point.y][f_point.z-1];
            int wire_length = 1;
            auto new_pos = new_gcell.passingNets.find(NetId);
            if(new_pos != new_gcell.passingNets.end() && new_pos->second != 0)
                wire_length = 0;
            int remain = new_gcell.capacity-new_gcell.demand-wire_length;
            if(remain > 0 && visited_p.find(new_p) == visited_p.end()) {
                p_q.emplace(new_p, new_gcell.demand+wire_length+distance(new_p,sink));
                visited_p[new_p] = f_point;
            }
        }
    }
    return find_flag;
}

int RoutingGraph::tree2tree_routing(priority_queue<pair<Point,int>, vector<pair<Point,int>>, greater<pair<Point,int>>>& p_q, Point b_min, 
    Point b_max, unordered_set<int>& source_comp_set, unordered_set<int>& sink_comp_set, 
    vector<vector<vector<int>>>& comp_grid_map, vector<vector<vector<int>>>& cost_grid_map, int NetId, 
    unordered_map<Point,Point,MyHashFunction>& visited_p, Point& reach_p) {
    while(!p_q.empty()) {
        auto frontier = p_q.top();
        p_q.pop();
        auto& f_point = frontier.first;
        int cost = frontier.second;
        Point local_p = f_point-Point(b_min.x,b_min.y,b_min.z);
        if(local_p.x<0 || local_p.y<0 || local_p.z<0 || f_point.x>b_max.x || f_point.y>b_max.y || f_point.z>b_max.z)
            continue;
        int tmp_comp = comp_grid_map[local_p.x][local_p.y][local_p.z];
        int via_cost = 0;
        //cout << "f_point: " << f_point << " " << tmp_comp << endl;
        if(sink_comp_set.find(tmp_comp) != sink_comp_set.end()){
            // find sink
            reach_p = f_point;
            return (nets[NetId].branch_nodes.find(reach_p) != nets[NetId].branch_nodes.end()) ? 1 : 2;
        }
        // up
        if(f_point.z%2==1 && f_point.x < b_max.x) {
            Point new_p(f_point.x+1,f_point.y,f_point.z);
            Gcell& new_gcell = grids[f_point.x+1][f_point.y][f_point.z];
            int wire_length = 1;
            auto new_pos = new_gcell.passingNets.find(NetId);
            if(new_pos != new_gcell.passingNets.end() && new_pos->second != 0)
                wire_length = 0;
            int remain = new_gcell.capacity-new_gcell.demand-wire_length;
            int& past_cost = cost_grid_map[local_p.x+1][local_p.y][local_p.z];
            //int future_cost = cost + 1;
            int future_cost;
            if(case3) future_cost = cost + 1;
            else future_cost = (remain > new_gcell.capacity/6) ? cost + 1 : cost + 1 + 1;
            //int future_cost = (remain > new_gcell.capacity/6) ? cost + 1 : cost + 1 + 1;
            if(remain > 0 && visited_p.find(new_p) == visited_p.end() && past_cost > future_cost) {    
                p_q.emplace(new_p, future_cost);
                past_cost = future_cost;
                visited_p[new_p] = f_point;
            }
        }
        // down
        if(f_point.z%2==1 && f_point.x > b_min.x) {
            Point new_p(f_point.x-1,f_point.y,f_point.z);
            Gcell& new_gcell = grids[f_point.x-1][f_point.y][f_point.z];
            int wire_length = 1;
            auto new_pos = new_gcell.passingNets.find(NetId);
            if(new_pos != new_gcell.passingNets.end() && new_pos->second != 0)
                wire_length = 0;
            int remain = new_gcell.capacity-new_gcell.demand-wire_length;
            int& past_cost = cost_grid_map[local_p.x-1][local_p.y][local_p.z];
            //int future_cost = cost + 1;
            int future_cost;
            if(case3) future_cost = cost + 1;
            else future_cost = (remain > new_gcell.capacity/6) ? cost + 1 : cost + 1 + 1;
            //int future_cost = (remain > new_gcell.capacity/6) ? cost + 1 : cost + 1 + 1;
            if(remain > 0 && visited_p.find(new_p) == visited_p.end() && past_cost > future_cost) {       
                p_q.emplace(new_p, future_cost);
                past_cost = future_cost;
                visited_p[new_p] = f_point;
            }
        }
        // right
        if(f_point.z%2==0 && f_point.y < b_max.y) {
            Point new_p(f_point.x,f_point.y+1,f_point.z);
            Gcell& new_gcell = grids[f_point.x][f_point.y+1][f_point.z];
            int wire_length = 1;
            auto new_pos = new_gcell.passingNets.find(NetId);
            if(new_pos != new_gcell.passingNets.end() && new_pos->second != 0)
                wire_length = 0;
            int remain = new_gcell.capacity-new_gcell.demand-wire_length;
            int& past_cost = cost_grid_map[local_p.x][local_p.y+1][local_p.z];
            //int future_cost = cost + 1;
            int future_cost;
            if(case3) future_cost = cost + 1;
            else future_cost = (remain > new_gcell.capacity/6) ? cost + 1 : cost + 1 + 1;
            //int future_cost = (remain > new_gcell.capacity/6) ? cost + 1 : cost + 1 + 1;
            if(remain > 0 && visited_p.find(new_p) == visited_p.end() && past_cost > future_cost) {       
                p_q.emplace(new_p, future_cost);
                past_cost = future_cost;
                visited_p[new_p] = f_point;
            }
        }
        // left
        if(f_point.z%2==0 && f_point.y > b_min.y) {
            Point new_p(f_point.x,f_point.y-1,f_point.z);
            Gcell& new_gcell = grids[f_point.x][f_point.y-1][f_point.z];
            int wire_length = 1;
            auto new_pos = new_gcell.passingNets.find(NetId);
            if(new_pos != new_gcell.passingNets.end() && new_pos->second != 0)
                wire_length = 0;
            int remain = new_gcell.capacity-new_gcell.demand-wire_length;
            int& past_cost = cost_grid_map[local_p.x][local_p.y-1][local_p.z];
            //int future_cost = cost + 1;
            int future_cost;
            if(case3) future_cost = cost + 1;
            else future_cost = (remain > new_gcell.capacity/6) ? cost + 1 : cost + 1 + 1;
            //int future_cost = (remain > new_gcell.capacity/6) ? cost + 1 : cost + 1 + 1;
            if(remain > 0 && visited_p.find(new_p) == visited_p.end() && past_cost > future_cost) {       
                p_q.emplace(new_p, future_cost);
                past_cost = future_cost;
                visited_p[new_p] = f_point;
            }
        }
        // top
        if(f_point.z < b_max.z) {
            Point new_p(f_point.x,f_point.y,f_point.z+1);
            Gcell& new_gcell = grids[f_point.x][f_point.y][f_point.z+1];
            int wire_length = 1;
            auto new_pos = new_gcell.passingNets.find(NetId);
            if(new_pos != new_gcell.passingNets.end() && new_pos->second != 0)
                wire_length = 0;
            int remain = new_gcell.capacity-new_gcell.demand-wire_length;
            int& past_cost = cost_grid_map[local_p.x][local_p.y][local_p.z+1];
            //int future_cost = cost + 1;
            int future_cost;
            if(case3) future_cost = cost + 1;
            else future_cost = (remain > new_gcell.capacity/6) ? cost + 1 : cost + 1 + 1;
            //int future_cost = (remain > new_gcell.capacity/6) ? cost + 1 : cost + 1 + 1;
            if(remain > 0 && visited_p.find(new_p) == visited_p.end() && past_cost > future_cost + via_cost) {       
                p_q.emplace(new_p, future_cost + via_cost);
                past_cost = future_cost+via_cost;
                visited_p[new_p] = f_point;
            }
        }
        // down
        if(f_point.z > b_min.z) {
            Point new_p(f_point.x,f_point.y,f_point.z-1);
            Gcell& new_gcell = grids[f_point.x][f_point.y][f_point.z-1];
            int wire_length = 1;
            auto new_pos = new_gcell.passingNets.find(NetId);
            if(new_pos != new_gcell.passingNets.end() && new_pos->second != 0)
                wire_length = 0;
            int remain = new_gcell.capacity-new_gcell.demand-wire_length;
            int& past_cost = cost_grid_map[local_p.x][local_p.y][local_p.z-1];
            //int future_cost = cost + 1;
            //int future_cost;
            //if(case3) future_cost = cost + 1;
            //else future_cost = (remain > new_gcell.capacity/6) ? cost + 1 : cost + 1 + 1;
            int future_cost = (remain > new_gcell.capacity/6) ? cost + 1 : cost + 1 + 1;
            if(remain > 0 && visited_p.find(new_p) == visited_p.end() && past_cost > future_cost + via_cost) {       
                p_q.emplace(new_p, future_cost + via_cost);
                past_cost = future_cost+via_cost;
                visited_p[new_p] = f_point;
            }
        }
    }
    return 0;
}

void RoutingGraph::add_component_in_pq(priority_queue<pair<Point,int>, vector<pair<Point,int>>, greater<pair<Point,int>>>& p_q, int source_comp, 
    unordered_map<Point, int, MyHashFunction>& component_map, int netId, unordered_map<Point,Point,MyHashFunction>& visited_p) {
    unordered_set<Point, MyHashFunction> sink_set;
    auto& net = nets[netId];
    for(auto& comp : component_map) {
        if(comp.second == source_comp) {
            // handle point there don't have neighbor, also should add into p_q
            if(net.branch_nodes[comp.first].neighbors.size() == 0) {
                Gcell& gcell = grids[comp.first.x][comp.first.y][comp.first.z];           
                auto pos = gcell.passingNets.find(netId);
                int wire_length = (pos != gcell.passingNets.end() && pos->second != 0) ? 0 : 1;
                int remain = gcell.capacity-gcell.demand-wire_length;
                if(remain >= 0) {
                    p_q.emplace(comp.first, 0);
                }
                return;
            }   
            for(auto& neighbor : net.branch_nodes[comp.first].neighbors) {
                if(sink_set.find(neighbor.first) != sink_set.end())
                    continue;
                // push path comp tag in grid
                for(auto& path : neighbor.second.paths)
                    add_path_comp_in_pq(p_q, path, netId, visited_p);
                sink_set.insert(comp.first);
            }
        }
    }
    
}

bool RoutingGraph::connect_all_nets(unordered_map<int, vector<tuple<Point,Point,int,Node,TwoPinNet>> >& open_nets, int& net_wirelength, 
    unordered_map<int, int>& netK, vector<pair<Point, int>>& point_nets) {
    for(auto& net_open_net : open_nets) {
        auto& net = nets[net_open_net.first];
        //cout << "net " << net.netId << endl;
        unordered_map<Point,Point,MyHashFunction> visited_p;
        unordered_map<Point, int, MyHashFunction> component_map;
        net.set_point_component(component_map);
        unordered_set<int> source_comp_set, sink_comp_set; 
        // set sink_comp_set
        for(auto& comp : component_map)
            sink_comp_set.insert(comp.second);
        Point b_min(INT32_MAX, INT32_MAX, INT32_MAX), b_max(0,0,0);
        net.calc_bounding_box(b_min, b_max, this->layer);
        //cout << "b_box " << b_min << " " << b_max << endl;
        if(b_min.z > b_max.z)  return 0;
        vector<vector<vector<int>>> comp_grid_map, cost_grid_map;
        set_all_comp_grid_map(comp_grid_map, cost_grid_map, net.netId, component_map, b_min, b_max);
        int source_comp = component_map.begin()->second;
        priority_queue<pair<Point,int>, vector<pair<Point,int>>, greater<pair<Point,int>>> p_q;
        while(sink_comp_set.size() > 1) {
            add_component_in_pq(p_q, source_comp, component_map, net.netId, visited_p);
            source_comp_set.insert(source_comp);
            sink_comp_set.erase(source_comp);
            Point reach_p, source;
            int find_flg = tree2tree_routing(p_q, b_min, b_max, source_comp_set, sink_comp_set, comp_grid_map, cost_grid_map, net.netId, visited_p, reach_p);         
            if(find_flg == 0) {
                // reverse
                return 0;
            }
            // success one net
            source_comp = comp_grid_map[reach_p.x-b_min.x][reach_p.y-b_min.y][reach_p.z-b_min.z];   
            //cout << "#find flg: " << find_flg << " reach_p: " << reach_p << " comp " << source_comp << endl;
            TwoPinNet two_pin = convert_path_to_twopin_t2t(source, reach_p, b_min, visited_p, source_comp_set, comp_grid_map);
            if(source == reach_p) {
                net.branch_nodes[reach_p].node.type = 1;
                continue;
            }
            net_wirelength -= two_pin.wire_length;
            if(net_wirelength <= 0) {
                // worse routing
                return 0;
            }
            if(find_flg == 2) {
                if(net.branch_nodes.find(reach_p) == net.branch_nodes.end())
                    net.insert_steiner_point(reach_p, grids);
            }
            if(net.branch_nodes.find(source) == net.branch_nodes.end())
                net.insert_steiner_point(source, grids);
            point_nets.emplace_back(reach_p, net.netId);
            // rebuild branch_nodes
            net.branch_nodes[source].neighbors.emplace_back(reach_p,two_pin);
            net.branch_nodes[reach_p].neighbors.emplace_back(source,two_pin_reverse(two_pin));
            net.wire_length += two_pin.wire_length;
            // add two_pin demand into graph
            net.add_twopin_demand_into_graph(two_pin, grids);
            for(auto& path : two_pin.paths)
                add_path_comp_in_comp_grid(comp_grid_map, path, b_min, source_comp);       
            netK[net.netId]++;
        }
    }
    return 1;
}

void Net::calc_bounding_box(Point& min, Point& max, int max_layer) {
    for(auto& node : branch_nodes) {
        min.x = (min.x < node.first.x) ? min.x : node.first.x;
        min.y = (min.y < node.first.y) ? min.y : node.first.y;
        min.z = (min.z < node.first.z) ? min.z : node.first.z;
        max.x = (max.x > node.first.x) ? max.x : node.first.x;
        max.y = (max.y > node.first.y) ? max.y : node.first.y;
        max.z = (max.z > node.first.z) ? max.z : node.first.z;
    }
    max.z = (max.z<max_layer-1) ? max.z+1 : max.z;
}

int RoutingGraph::A_star_pin2component_routing(Point source, Point sink, int NetId, unordered_map<Point,Point,MyHashFunction>& visited_p,
    unordered_map<Point, int, MyHashFunction>& component_map, Point& reach_p) {
    int min_x = min(source.x, sink.x);
    int max_x = max(source.x, sink.x);
    int min_y = min(source.y, sink.y);
    int max_y = max(source.y, sink.y);
    int min_l = min(source.z, sink.z);
    int max_l = max(source.z, sink.z);
    visited_p[source] = source;
    if(source == sink) return 0;
    max_l = (max_l<this->layer-1) ? max_l+1 : max_l;
    int min_r_l = nets[NetId].minRoutingLayer;
    if(min_r_l>0) return 0;
    Point ori_source = source, ori_sink = sink;
    if(min_l < min_r_l) {
        min_l = min_r_l;
        if(source.z < min_r_l) {
            source.z = min_r_l;
            visited_p[source] = ori_source;
        }
        if(sink.z < min_r_l) {
            sink.z = min_r_l;
            visited_p[ori_sink] = sink;
        }
    }
    if(max_l < min_l) return 0;
    priority_queue<pair<Point,int>, vector<pair<Point,int>>, greater<pair<Point,int>>> p_q;
    // initial comp_grid_map
    vector<vector<vector<int>>> comp_grid_map;
    set_comp_grid_map(comp_grid_map, NetId, sink, component_map, Point(min_x,min_y,min_l), Point(max_x,max_y,max_l));

    Gcell& gcell = grids[source.x][source.y][source.z];
    int wire_length = 1;
    auto pos = gcell.passingNets.find(NetId);
    if(pos != gcell.passingNets.end() && pos->second != 0)
        wire_length = 0;
    int remain = gcell.capacity-gcell.demand-wire_length;
    if(remain < 0)
        return 0;
    p_q.emplace(source, gcell.demand+wire_length+distance(source,sink));
    while(!p_q.empty()) {
        auto frontier = p_q.top();
        p_q.pop();
        auto& f_point = frontier.first;
        if(f_point == sink) {
            reach_p = sink;
            return 1;
        }
        Point local_p = f_point-Point(min_x,min_y,min_l);
        if(comp_grid_map[local_p.x][local_p.y][local_p.z] > 0) {
            reach_p = f_point;
            return 2;
        }
        // up
        if(f_point.z%2==1 && f_point.x < max_x) {
            Point new_p(f_point.x+1,f_point.y,f_point.z);
            Gcell& new_gcell = grids[f_point.x+1][f_point.y][f_point.z];
            int wire_length = 1;
            auto new_pos = new_gcell.passingNets.find(NetId);
            if(new_pos != new_gcell.passingNets.end() && new_pos->second != 0)
                wire_length = 0;
            int remain = new_gcell.capacity-new_gcell.demand-wire_length;
            if(remain > 0 && visited_p.find(new_p) == visited_p.end()) {              
                p_q.emplace(new_p, new_gcell.demand+wire_length+distance(new_p,sink));
                visited_p[new_p] = f_point;
            }
        }
        // down
        if(f_point.z%2==1 && f_point.x > min_x) {
            Point new_p(f_point.x-1,f_point.y,f_point.z);
            Gcell& new_gcell = grids[f_point.x-1][f_point.y][f_point.z];
            int wire_length = 1;
            auto new_pos = new_gcell.passingNets.find(NetId);
            if(new_pos != new_gcell.passingNets.end() && new_pos->second != 0)
                wire_length = 0;
            int remain = new_gcell.capacity-new_gcell.demand-wire_length;
            if(remain > 0 && visited_p.find(new_p) == visited_p.end()) {
                p_q.emplace(new_p, new_gcell.demand+wire_length+distance(new_p,sink));
                visited_p[new_p] = f_point;
            }
        }
        // right
        if(f_point.z%2==0 && f_point.y < max_y) {
            Point new_p(f_point.x,f_point.y+1,f_point.z);
            Gcell& new_gcell = grids[f_point.x][f_point.y+1][f_point.z];
            int wire_length = 1;
            auto new_pos = new_gcell.passingNets.find(NetId);
            if(new_pos != new_gcell.passingNets.end() && new_pos->second != 0)
                wire_length = 0;
            int remain = new_gcell.capacity-new_gcell.demand-wire_length;
            if(remain > 0 && visited_p.find(new_p) == visited_p.end()) {
                p_q.emplace(new_p, new_gcell.demand+wire_length+distance(new_p,sink));
                visited_p[new_p] = f_point;
            }
        }
        // left
        if(f_point.z%2==0 && f_point.y > min_y) {
            Point new_p(f_point.x,f_point.y-1,f_point.z);
            Gcell& new_gcell = grids[f_point.x][f_point.y-1][f_point.z];
            int wire_length = 1;
            auto new_pos = new_gcell.passingNets.find(NetId);
            if(new_pos != new_gcell.passingNets.end() && new_pos->second != 0)
                wire_length = 0;
            int remain = new_gcell.capacity-new_gcell.demand-wire_length;
            if(remain > 0 && visited_p.find(new_p) == visited_p.end()) {
                p_q.emplace(new_p, new_gcell.demand+wire_length+distance(new_p,sink));
                visited_p[new_p] = f_point;
            }
        }
        // top
        if(f_point.z < max_l) {
            Point new_p(f_point.x,f_point.y,f_point.z+1);
            Gcell& new_gcell = grids[f_point.x][f_point.y][f_point.z+1];
            int wire_length = 1;
            auto new_pos = new_gcell.passingNets.find(NetId);
            if(new_pos != new_gcell.passingNets.end() && new_pos->second != 0)
                wire_length = 0;
            int remain = new_gcell.capacity-new_gcell.demand-wire_length;
            if(remain > 0 && visited_p.find(new_p) == visited_p.end()) {
                p_q.emplace(new_p, new_gcell.demand+wire_length+distance(new_p,sink));
                visited_p[new_p] = f_point;
            }
        }
        // down
        if(f_point.z > min_l) {
            Point new_p(f_point.x,f_point.y,f_point.z-1);
            Gcell& new_gcell = grids[f_point.x][f_point.y][f_point.z-1];
            int wire_length = 1;
            auto new_pos = new_gcell.passingNets.find(NetId);
            if(new_pos != new_gcell.passingNets.end() && new_pos->second != 0)
                wire_length = 0;
            int remain = new_gcell.capacity-new_gcell.demand-wire_length;
            if(remain > 0 && visited_p.find(new_p) == visited_p.end()) {
                p_q.emplace(new_p, new_gcell.demand+wire_length+distance(new_p,sink));
                visited_p[new_p] = f_point;
            }
        }
    }
    return 0;
}

void RoutingGraph::set_comp_grid_map(std::vector<std::vector<std::vector<int>>>& comp_grid_map, int netId, Point sink,
    unordered_map<Point, int, MyHashFunction>& component_map, Point box_min, Point box_max) {
    comp_grid_map.resize(box_max.x-box_min.x+1);
    for(auto& comp_x : comp_grid_map) {
        comp_x.resize(box_max.y-box_min.y+1);
        for(auto& comp_y : comp_x)
            comp_y.resize(box_max.z-box_min.z+1, 0);
    }
    Net& net = nets[netId];
    int sink_comp = component_map[sink];
    unordered_set<Point, MyHashFunction> sink_set;
    for(auto& comp : component_map) {
        if(comp.second == sink_comp) {
            for(auto& neighbor : net.branch_nodes[comp.first].neighbors) {
                if(sink_set.find(neighbor.first) != sink_set.end())
                    continue;
                // push path comp tag in grid
                for(auto& path : neighbor.second.paths)
                    add_path_comp_in_comp_grid(comp_grid_map, path, box_min, sink_comp);
                sink_set.insert(comp.first);
            }
        }
    }
}

void RoutingGraph::set_all_comp_grid_map(std::vector<std::vector<std::vector<int>>>& comp_grid_map, 
    std::vector<std::vector<std::vector<int>>>& cost_grid_map, int netId,
    unordered_map<Point, int, MyHashFunction>& component_map, Point box_min, Point box_max) {
    comp_grid_map.resize(box_max.x-box_min.x+1);
    cost_grid_map.resize(box_max.x-box_min.x+1);
    for(int x=0; x<comp_grid_map.size(); x++) {
        auto& comp_x = comp_grid_map[x];
        auto& cost_x = cost_grid_map[x];
        comp_x.resize(box_max.y-box_min.y+1);
        cost_x.resize(box_max.y-box_min.y+1);
        for(int y=0; y<comp_x.size(); y++) {
            auto& comp_y = comp_x[y];
            auto& cost_y = cost_x[y];
            comp_y.resize(box_max.z-box_min.z+1, 0);
            cost_y.resize(box_max.z-box_min.z+1, INT32_MAX);
        }         
    }
    Net& net = nets[netId];
    unordered_set<Point, MyHashFunction> treversed;
    for(auto& comp : component_map) {
        // handle node wi 
        Point local_p = comp.first - box_min;
        if(local_p.x<0 || local_p.y<0 || local_p.z<0) continue;
        if(local_p.x >= comp_grid_map.size() || local_p.y >= comp_grid_map[0].size() || local_p.z >= comp_grid_map[0][0].size())
            continue;
        comp_grid_map[local_p.x][local_p.y][local_p.z] = comp.second;

        for(auto& neighbor : net.branch_nodes[comp.first].neighbors) {
            if(treversed.find(neighbor.first) != treversed.end())
                continue;
            // push path comp tag in grid
            for(auto& path : neighbor.second.paths)
                add_path_comp_in_comp_grid(comp_grid_map, path, box_min, comp.second);             
            treversed.insert(comp.first);
        }
    }
}

void RoutingGraph::add_path_comp_in_comp_grid(std::vector<std::vector<std::vector<int>>>& comp_grid_map, 
    pair<Point, Point> path, Point box_min, int sink_comp) {
    if( !(path.first<=path.second) )
        swap(path.first, path.second);
    Point dir = norm(path.second-path.first);
    for(Point p = path.first; p <= path.second; p = p+dir) {
        Point local_p = p - box_min;
        if(local_p.x<0 || local_p.y<0 || local_p.z<0) continue;
        if(local_p.x >= comp_grid_map.size() || local_p.y >= comp_grid_map[0].size() || local_p.z >= comp_grid_map[0][0].size())
            continue;
        comp_grid_map[local_p.x][local_p.y][local_p.z] = sink_comp;
    }
}

void RoutingGraph::add_path_comp_in_pq(priority_queue<pair<Point,int>, vector<pair<Point,int>>, greater<pair<Point,int>>>& p_q, pair<Point, Point> path, int netId, 
    unordered_map<Point,Point,MyHashFunction>& visited_p) {
    if( !(path.first<=path.second) )
        swap(path.first, path.second);
    Point dir = norm(path.second-path.first);
    for(Point p = path.first; p <= path.second; p = p+dir) {
        Gcell& gcell = grids[p.x][p.y][p.z];
        int wire_length = 1;
        auto pos = gcell.passingNets.find(netId);
        if(pos != gcell.passingNets.end() && pos->second != 0)
            wire_length = 0;
        int remain = gcell.capacity-gcell.demand-wire_length;
        if(remain < 0)
            continue;
        p_q.emplace(p, 0);
        //cout << "p_q add: " << p << endl;
        visited_p[p] = p;
    }
}

void RoutingGraph::reroute_all_net() {
    for(auto it = nets.begin(); it != nets.end(); ++it) {
        Net& net = *it;
        unordered_map<Point, bool, MyHashFunction> visited;
        if(net.minRoutingLayer > 0) continue;
        for(auto _it = net.branch_nodes.begin(); _it != net.branch_nodes.end(); ++_it) {
            const Point& start = _it->first;
            TreeNode& treeNode = _it->second;
            visited[start] = true;
            for(auto __it = treeNode.neighbors.begin(); __it != treeNode.neighbors.end(); ++__it) {
                Point& end = __it->first;
                if(visited[end] == false) {
                    net.del_twoPinNet_from_graph(__it->second, grids);
                    unordered_map<Point,Point,MyHashFunction> visited_p;
                    A_star_routing(start, end, net.netId, visited_p);
                    Point tmp_node = end;
                    int dir = 0; // 1: vertical, 2: horizontal, 3:via
                    TwoPinNet two_pin;
                    two_pin.n1.p = start;
                    two_pin.n2.p = end;
                    two_pin.paths.push_back({end,visited_p[end]});
                    // build TwopinNet
                    while(true) {
                        Point back_p = visited_p[tmp_node];
                        int new_dir = 0;
                        if(back_p.x != tmp_node.x) new_dir = 1;
                        if(back_p.y != tmp_node.y) new_dir = 2;
                        if(back_p.z != tmp_node.z) new_dir = 3;
                        if(new_dir != dir) {
                            two_pin.paths.back().second = tmp_node;
                            if(tmp_node!=start && tmp_node != end)
                                two_pin.paths.push_back({tmp_node,back_p});
                            dir = new_dir;
                        }
                        if(tmp_node == start)
                            break;
                        tmp_node = back_p;
                    }
                    // rebuild branch_nodes
                    for(int i = 0; i < net.branch_nodes[start].neighbors.size(); ++i) {
                        if(net.branch_nodes[start].neighbors[i].first == end) {
                            net.branch_nodes[start].neighbors[i].second = two_pin;
                            break;
                        }
                    }
                    for(int i = 0; i < net.branch_nodes[end].neighbors.size(); ++i) {
                        if(net.branch_nodes[end].neighbors[i].first == start) {
                            net.branch_nodes[end].neighbors[i].second = two_pin;
                            break;
                        }
                    }
                    // add two_pin demand into graph
                }
            }
        }
    }
}

// void RoutingGraph::reroute_cell_two_pin_net(int net_Id) {
//     for(auto& node : nets[net_Id].branch_nodes) {
//         if(node.second.node.type == 1) {
            
//         }
//     }
// }

void TwoPinNet::update_wire_length() {
    wire_length = 0;
    for(auto& path : paths) {
        wire_length += abs((path.first.x-path.second.x)) + abs((path.first.y-path.second.y)) + abs((path.first.z-path.second.z));
    }
}

int RoutingGraph::get_edge_cost(Point& from, Point& to) {
    int alpha = 1, beta = 1;
    return alpha*max(grids[from.x][from.y][from.z].demand, grids[to.x][to.y][to.z].demand) + beta*1;
}

void RoutingGraph::swap_into_optimal_region(void) {
    vector<vector<pair<Point, int>>> cells_pos(cellInstances.size());
    map<pair<int,int>, bool> haveSwap;
    for(int cellIdx = 0; cellIdx < cellInstances.size(); cellIdx++) {
        Cell& cell = cellInstances[cellIdx];
        find_optimal_pos(cell, cells_pos[cellIdx]);
    }
    for(int cellIdx = 0; cellIdx < cellInstances.size(); cellIdx++) {
        if(movedCell.size() >= maxCellMove) return;
        Cell& cell = cellInstances[cellIdx];
        vector<pair<Point,int>>& cell_pos = cells_pos[cellIdx];
        for(auto it = cell_pos.begin(); it != cell_pos.end(); ++it) {
            Point& opt = it->first;
            //cannot use reference since swap will affect placement
            unordered_set<int> cellsInTheRegion = placement[opt.x][opt.y];
            for(int _cellIdx : cellsInTheRegion) {
                vector<pair<Point,int>>& candidateCellPos = cells_pos[_cellIdx];
                for(auto _it = candidateCellPos.begin(); _it != candidateCellPos.end(); ++_it) {
                    Point& candidateOpt = _it->first;
                    if(candidateOpt.x == cell.x && candidateOpt.y == cell.y && haveSwap[{cellIdx,_cellIdx}] == false) {
                        //cout << "Cell " << cellIdx << " and Cell " << _cellIdx << " can swap\n";
                        int wireLength = 0;
                        bool success = swap_two_cells(cellIdx, _cellIdx, wireLength);
                        if(success) {
                            haveSwap[{cellIdx,_cellIdx}] = true;
                            haveSwap[{_cellIdx,cellIdx}] = true;
                        }
                    }
                }
            }
        }
    }
}

bool RoutingGraph::swap_two_cells(int cell_idx1, int cell_idx2, int& net_wirelength) {
    //cellIdx -> netId -> number to delete
    unordered_map<int, unordered_map<int, int>> netK;
    Cell& cell1 = cellInstances[cell_idx1];
    Cell& cell2 = cellInstances[cell_idx2];
    //printf("Cell1: %d (%d,%d) Cell2: %d (%d,%d)\n", cell_idx1, cell1.x, cell1.y, cell_idx2, cell2.x, cell2.y);
    if(cell1.x == cell2.x && cell1.y == cell2.y) return 1;
    int cell1Exist = (movedCell.find(cell_idx1) == movedCell.end())? 1:0;
    int cell2Exist = (movedCell.find(cell_idx2) == movedCell.end())? 1:0;
    if((movedCell.size() + cell1Exist + cell2Exist) > maxCellMove) return 0;
    //if(check_cell_cost_in_graph(cell2.x, cell2.y, cell1) < 0 || check_cell_cost_in_graph(cell1.x, cell1.y, cell2) < 0) return 0;
    for(auto& pin : cell1.pins) {
        if(pin.connectedNet == -1) continue;
        Net& net = nets[pin.connectedNet];
        Point p(cell1.x, cell1.y, pin.layer);
        TreeNode& treeNode = net.branch_nodes[p];
        for(auto it = treeNode.neighbors.begin(); it != treeNode.neighbors.end(); ++it) {
            Point& neighbor = it->first;
            if(neighbor.x == cell2.x && neighbor.y == cell2.y) return 0;
        }
    }
    int cell1OriX = cell1.x, cell1OriY = cell1.y, cell2OriX = cell2.x, cell2OriY = cell2.y;
    // key is net, value is tuple
    // source, sink, netId, source.node, twopin
    unordered_map<int, vector<tuple<Point,Point,int,Node,TwoPinNet>>> open_nets1;
    unordered_map<int, vector<tuple<Point,Point,int,Node,TwoPinNet>>> open_nets2;
    // key is netId, value is branch_node
    unordered_map<int,unordered_map<Point, TreeNode, MyHashFunction>> branchs_copy;
    unordered_set<int> updated_branchs;
    unordered_map<int, map<tuple<int, int, int, int, int, int>, bool>> alreadyInOpenNets;
    int net_wirelength1 = 0, net_wirelength2 = 0;
    unordered_set<int> net_id_pool;
    //cout << "Cell: " << cell_idx << endl;
    for(auto& pin : cell1.pins) {
        if(pin.connectedNet == -1)
            continue;
        auto& net = nets[pin.connectedNet];
        net_id_pool.insert(net.netId);
        // copy branch_node for reverse
        if(branchs_copy.find(pin.connectedNet) == branchs_copy.end())
            branchs_copy[pin.connectedNet] = net.branch_nodes;
        Point cell_p(cell1.x, cell1.y, pin.layer);
        unordered_set<Point, MyHashFunction> used_local_node;
        for(auto& neighbor : net.branch_nodes[cell_p].neighbors) {
            Point neighbor_p = neighbor.first;
            // local pin can not put in open net twice
            if(used_local_node.find(neighbor_p) != used_local_node.end()) {
                cout << "skip\n";
                continue;
            }
            if(net.branch_nodes[neighbor_p].node.mergedLocalPins.size() >= 2)
                used_local_node.insert(neighbor_p);
            if(alreadyInOpenNets[net.netId][{cell2.x,cell2.y,pin.layer,neighbor_p.x, neighbor_p.y, neighbor_p.z}] == true) continue;
            neighbor.second.update_wire_length();
            net_wirelength1 += neighbor.second.wire_length;
            open_nets1[net.netId].emplace_back(Point(cell2.x,cell2.y,pin.layer), neighbor_p, net.netId, net.branch_nodes[cell_p].node, neighbor.second);
            alreadyInOpenNets[net.netId][{cell2.x,cell2.y,pin.layer,neighbor_p.x, neighbor_p.y, neighbor_p.z}] = true;
        }
        if(net.branch_nodes[cell_p].node.type == 2 && net.branch_nodes[cell_p].neighbors.empty()
        && alreadyInOpenNets[net.netId][{cell2.x,cell2.y,pin.layer,cell1OriX, cell1OriY, pin.layer}] == false) {
            net_wirelength1++;
            open_nets1[net.netId].emplace_back(Point(cell2.x,cell2.y,pin.layer), Point(cell1OriX,cell1OriY,pin.layer), net.netId, Node(), TwoPinNet());
            alreadyInOpenNets[net.netId][{cell2.x,cell2.y,pin.layer,cell1OriX, cell1OriY, pin.layer}] = true;
        }
    }
    for(auto& pin : cell2.pins) {
        if(pin.connectedNet == -1) continue;
        Net& net = nets[pin.connectedNet];
        net_id_pool.insert(net.netId);
        if(branchs_copy.find(pin.connectedNet) == branchs_copy.end()) branchs_copy[pin.connectedNet] = net.branch_nodes;
        Point cell_p(cell2.x, cell2.y, pin.layer);
        unordered_set<Point, MyHashFunction> used_local_node;
        for(auto& neighbor : net.branch_nodes[cell_p].neighbors) {
            Point neighbor_p = neighbor.first;
            if(used_local_node.find(neighbor_p) != used_local_node.end()) {
                cout << "skip\n";
                continue;
            }
            if(net.branch_nodes[neighbor_p].node.mergedLocalPins.size() >= 2) used_local_node.insert(neighbor_p);
            if(alreadyInOpenNets[net.netId][{cell1.x,cell1.y,pin.layer,neighbor_p.x, neighbor_p.y, neighbor_p.z}] == true) continue;
            neighbor.second.update_wire_length();
            net_wirelength2 += neighbor.second.wire_length;
            open_nets2[net.netId].emplace_back(Point(cell1.x,cell1.y,pin.layer), neighbor_p, net.netId, net.branch_nodes[cell_p].node, neighbor.second);
            alreadyInOpenNets[net.netId][{cell1.x,cell1.y,pin.layer,neighbor_p.x, neighbor_p.y, neighbor_p.z}] = true;
        }
        if(net.branch_nodes[cell_p].node.type == 2 && net.branch_nodes[cell_p].neighbors.empty()
        && alreadyInOpenNets[net.netId][{cell1.x,cell1.y,pin.layer,cell2OriX, cell2OriY, pin.layer}] == false) {
            net_wirelength2++;
            open_nets2[net.netId].emplace_back(Point(cell1.x,cell1.y,pin.layer), Point(cell2OriX,cell2OriY,pin.layer), net.netId, Node(), TwoPinNet());
            alreadyInOpenNets[net.netId][{cell1.x,cell1.y,pin.layer,cell2OriX, cell2OriY, pin.layer}] = true;
        }
    }
    int all_net_wl = 0;
    for(auto& netId : net_id_pool) {
        nets[netId].update_wirelength();
        all_net_wl += nets[netId].wire_length;
    }
    del_cell_neighbor(cell_idx1);
    del_cell_neighbor(cell_idx2);
    add_cell(cell2OriX, cell2OriY, cell_idx1);
    add_cell(cell1OriX, cell1OriY, cell_idx2);
    vector<pair<Point, int>> point_nets1, point_nets2;
    bool routing_success1 = connect_all_nets(open_nets1, net_wirelength1, netK[cell_idx1], point_nets1);
    bool routing_success2 = connect_all_nets(open_nets2, net_wirelength2, netK[cell_idx2], point_nets2);
    if(routing_success1 && routing_success2) {
        for(auto& netId : net_id_pool) {
            nets[netId].update_wirelength();
            all_net_wl -= nets[netId].wire_length;
        }
        if(all_net_wl <= 0) routing_success1 = 0;
        else net_wirelength = all_net_wl;
        for(int i = 0; i < layer; ++i) {
            if(grids[cell1OriX][cell1OriY][i].get_remaining() < 0) {
                routing_success1 = 0;
                break;
            }
        }
        for(int i = 0; i < layer; ++i) {
            if(grids[cell2OriX][cell2OriY][i].get_remaining() < 0) {
                routing_success1 = 0;
                break;
            }
        }
        if(cell1OriY > 0) {
            for(int i = 0; i < layer; ++i) {
                if(grids[cell1OriX][cell1OriY-1][i].get_remaining() < 0) {
                    routing_success1 = 0;
                    break;
                }
            }
        }
        if(cell1OriY < (column-1)) {
            for(int i = 0; i < layer; ++i) {
                if(grids[cell1OriX][cell1OriY+1][i].get_remaining() < 0) {
                    routing_success1 = 0;
                    break;
                }
            }
        }
        if(cell2OriY > 0) {
            for(int i = 0; i < layer; ++i) {
                if(grids[cell2OriX][cell2OriY-1][i].get_remaining() < 0) {
                    routing_success1 = 0;
                    break;
                }
            }
        }
        if(cell2OriY < (column-1)) {
            for(int i = 0; i < layer; ++i) {
                if(grids[cell2OriX][cell2OriY+1][i].get_remaining() < 0) {
                    routing_success1 = 0;
                    break;
                }
            }
        }
    }
    if(routing_success1 && routing_success2) {
        movedCell.insert(cell_idx1);
        movedCell.insert(cell_idx2);
        return 1;
    }
    else {
        net_wirelength = 0;
        // reverse original net
        //del_cell_last_k_neighbor(cell_idx, netK);
        del_cell_point_net(cell_idx1, point_nets1);
        del_cell_point_net(cell_idx2, point_nets2);
        if(cell1OriX == cell1.originalX && cell1OriY == cell1.originalY)
            movedCell.erase(cell_idx1);
        if(cell2OriX == cell2.originalX && cell2OriY == cell2.originalY)
            movedCell.erase(cell_idx2);
        add_cell(cell1OriX, cell1OriY, cell_idx1);
        add_cell(cell2OriX, cell2OriY, cell_idx2);
        for(auto& net_open_net : open_nets1) {
            for(auto& open_net : net_open_net.second) {
                auto& net = nets[get<2>(open_net)];
                // resume net.branch_nodes
                if(updated_branchs.find(net.netId) == updated_branchs.end()) {
                    updated_branchs.insert(net.netId);
                    net.branch_nodes = branchs_copy[net.netId];
                }
                Point p1 = Point(cell1OriX, cell1OriY, get<0>(open_net).z);
                if(p1 == get<1>(open_net))
                    continue;
                // add two_pin demand into graph
                net.add_twopin_demand_into_graph(get<4>(open_net), grids);
            }
        }
        for(auto& net_open_net : open_nets2) {
            for(auto& open_net : net_open_net.second) {
                auto& net = nets[get<2>(open_net)];
                // resume net.branch_nodes
                if(updated_branchs.find(net.netId) == updated_branchs.end()) {
                    updated_branchs.insert(net.netId);
                    net.branch_nodes = branchs_copy[net.netId];
                }
                Point p1 = Point(cell2OriX, cell2OriY, get<0>(open_net).z);
                if(p1 == get<1>(open_net))
                    continue;
                // add two_pin demand into graph
                net.add_twopin_demand_into_graph(get<4>(open_net), grids);
            }
        }
        //cout << "#after end\n";
        return 0;
    }
}
