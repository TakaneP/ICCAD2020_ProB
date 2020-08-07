#include<bits/stdc++.h>
#include<unordered_set>

extern "C" {
#include"flute-3.1/flute.h"
}

#ifndef DATASTRUCTURE
#define DATASTRUCTURE

class SegmentTree;
struct Cell;
struct TreeNode;
struct Node;
struct TwoPinNet;

TwoPinNet two_pin_reverse(TwoPinNet two_pin);

struct Pin{
    Pin(): layer(-1), connectedNet(-1), pseudo(0) {}
    Pin(int l): layer(l), connectedNet(-1), pseudo(0) {}
    Pin(int l, bool p, int a): layer(l), pseudo(p), actualPinLayer(a) {}
    int layer;
    int connectedNet;
    bool pseudo;
    int actualPinLayer; //Only if pin is pseudom used to output
};

struct MasterCell{
    MasterCell() {}
    MasterCell(std::vector<Pin>& p, std::unordered_map<int, int>& b): pins(p), blockages(b) {}
    std::vector<Pin> pins;
    std::unordered_map<int, int> blockages; //layer's demand sum
};

struct Point{
    Point() {}
    Point(int i, int j, int k): x(i), y(j), z(k) {}
    friend std::ostream& operator<<(std::ostream& os, const Point& dt){
        os<<"("<<dt.x<<","<<dt.y<<","<<dt.z<<")";
        return os;
    }
    bool operator==(const Point& p2) const{
        return (this->x==p2.x) && (this->y==p2.y) && (this->z==p2.z);
    }
    bool operator!=(const Point& p2) const{
        return (this->x!=p2.x) || (this->y!=p2.y) || (this->z!=p2.z);
    }
    Point operator+(const Point& p2) const{
        return Point(this->x+p2.x, this->y+p2.y, this->z+p2.z);
    }
    Point operator-(const Point& p2) const{
        return Point(this->x-p2.x, this->y-p2.y, this->z-p2.z);
    }
    bool operator<=(const Point& p2) const{
        return (this->x<=p2.x) && (this->y<=p2.y) && (this->z<=p2.z);
    }
    int x,y,z;
};

class MyHashFunction { 
public: 
    size_t operator()(const Point& p) const;
    size_t operator()(const std::pair<Point,Point>& p) const;
}; 

struct Node{
    Point p;
    int type; //1: pin, 0 steiner node, -1 redundant point, 2: merged local pin
    std::vector<std::pair<int,int>> mergedLocalPins;
    bool operator==(const Node& p2){
        return this->p == p2.p;
    }
    Node() {
        p = Point(0,0,0);
        type = 0;
    }
    Node(Point point, int type): p(point), type(type) {}
};

struct TwoPinNet{
    Node n1, n2;
    std::vector<std::pair<Point,Point>> paths;
    int wire_length;
    void update_wire_length();
};

struct DegreeNode{
    // up/down is y axis, left/right is x axis, top/bottom is z axis
    DegreeNode() {up=0, down=0, left=0, right=0, top=0, bottom=0;}
    bool up, down, left, right, top, bottom;
    int return_degree(){
        return (int)up+down+left+right+top+bottom;
    }
};

struct TreeNode{
    TreeNode(){}
    TreeNode(Node n): node(n) {}
    Node node;
    std::vector<std::pair<Point,TwoPinNet>> neighbors;
};

struct Gcell;

struct Net{
    Net() {}
    int netId;
    int minRoutingLayer;
    std::vector<std::pair<int,int>> pins; //First: Cell Instance  Second: Pin
    std::vector<std::pair<Point,Point>> routingSegments;
    std::unordered_map <Point,TreeNode,MyHashFunction> branch_nodes;
    
    void convert_seg_to_2pin(std::vector<std::vector<std::vector<DegreeNode>>>& degreeMap, 
        std::vector<Cell>& cellInstances,
        std::vector<MasterCell>& masterCells,
        std::vector<std::vector<std::vector<Gcell>>>& grids
        );
    // add segment in passingMap, and construct steiner_map
    void set_passing_map(std::vector<std::vector<std::vector<DegreeNode>>>& degreeMap, std::vector<Cell>& cellInstances, 
        std::vector<MasterCell>& masterCells, std::unordered_set <Point,MyHashFunction>& pin_map, 
        std::unordered_set <Point,MyHashFunction>& steiner_map, int value);
    void traverse_passing_map(std::vector<std::vector<std::vector<DegreeNode>>>& degreeMap, 
        std::unordered_set <Point,MyHashFunction>& pin_map, 
        std::unordered_set <Point,MyHashFunction>& steiner_map,
        Point start_p,
        std::map<std::tuple<int,int,int>, std::vector<std::pair<int,int>>>& localNets,
        std::vector<TwoPinNet>& routingTree,
        std::vector<std::vector<std::vector<Gcell>>>& grids);
    Point return_next_dir(std::vector<std::vector<std::vector<DegreeNode>>>& degreeMap, Point now_p);
    bool check_map_legal(std::vector<std::vector<std::vector<DegreeNode>>>& degreeMap, Point now_p);
    bool check_map_dir(std::vector<std::vector<std::vector<DegreeNode>>>& degreeMap, Point now_p, Point dir);
    void decrese_degree_endpoint(std::vector<std::vector<std::vector<DegreeNode>>>& degreeMap, 
        Point now_p, Point dir);
    void decrese_degree_middle_p(std::vector<std::vector<std::vector<DegreeNode>>>& degreeMap, 
        Point now_p, Point dir);
    void print_two_pins(std::vector<TwoPinNet>& routingTree);
    void print_branch_nodes();
    void construct_branch_nodes(std::vector<TwoPinNet>& routingTree);
    void remove_dangling_wire(std::vector<std::vector<std::vector<Gcell>>>& grids);
    void clear_steiner_point(Point p, std::vector<std::vector<std::vector<Gcell>>>& grids);
    // construct MST to remove cycle in branch_nodes
    void remove_branch_cycle(std::vector<std::vector<std::vector<Gcell>>>& grids);
    void push_edge_in_queue(std::priority_queue<TwoPinNet, std::vector<TwoPinNet>, std::greater<TwoPinNet>>& frontier_edges);
    void del_net_from_graph(int x, int y, int z, std::vector<std::vector<std::vector<Gcell>>>& grids);
    void del_seg_demand(std::pair<Point,Point> segment, std::vector<std::vector<std::vector<Gcell>>>& grids);
    void del_seg_demand_from_graph(int x, int y, int z, std::vector<std::vector<std::vector<Gcell>>>& grids);
    void del_twoPinNet_from_graph(TwoPinNet& twoPinNet, std::vector<std::vector<std::vector<Gcell>>>& grids);
    void add_net_demand_into_graph(int x, int y, int z, std::vector<std::vector<std::vector<Gcell>>>& grids);
    void add_twopin_demand_into_graph(TwoPinNet& twoPinNet, std::vector<std::vector<std::vector<Gcell>>>& grids);
    bool insert_steiner_point(Point p);
    void set_point_component(std::unordered_map<Point, int, MyHashFunction>& component_map);
    void calc_bounding_box(Point& min, Point& max, int max_layer);
};

struct Cell : public MasterCell{
    Cell() {}
    Cell(std::vector<Pin>& p, std::unordered_map<int, int>& b, bool m, int _x, int _y, int mc): MasterCell(p,b), movable(m), x(_x), y(_y), mcType(mc) {
        originalX = x;
        originalY = y;
    }
    int mcType;
    bool movable;
    int x,y;
    int originalX, originalY;
    //std::unordered_set<int> connectedNets;    
};

struct Gcell{
    Gcell(): demand(0) {}
    Gcell(int c) : capacity(c), demand(0) {}     
    int get_remaining(void) {return capacity - demand;}
    int capacity;
    int demand;
    std::unordered_map<int, int> passingNets; //key: netId, value: sum of the number of pins and the number of segments in this Gcell for this net
};

class RoutingGraph{
public:
    RoutingGraph();
    ~RoutingGraph();
    void add_cell(int x, int y, int cellIndex);
    void del_cell(int cellIndex);
    void add_cell_demand_into_graph(int x, int y, int MCtype);
    void del_cell_demand_from_graph(int x, int y, int MCtype);
    void del_cell_neighbor(int cellIndex);
    void del_cell_last_k_neighbor(int cellIndex, std::unordered_map<int, int>& netK); //key: netId, value: delete last k neighbors for that net
    void construct_2pin_nets();
    void move_cells_force();
    void wirelength_driven_move();
    void swap_two_cells(int cell_idx1, int cell_idx2);
    bool move_cell_into_optimal_region(int cell_idx);
    void reroute_all_net();
    bool find_optimal_pos(Cell& cell, std::vector<std::pair<Point,int>>& cells_pos);
    // return cell profit after put in cell
    int check_cell_cost_in_graph(int x, int y, Cell& cell);
    int Z_shape_routing(Point source, Point sink, int NetId);
    void add_component_in_pq(std::priority_queue<std::pair<Point,int>>& p_q, int source_comp, 
        std::unordered_map<Point, int, MyHashFunction>& component_map, int netId, std::unordered_map<Point,Point,MyHashFunction>& visited_p);
    bool connect_all_nets(std::unordered_map<int, std::vector<std::tuple<Point,Point,int,Node,TwoPinNet>> >& open_nets, int& net_wirelength, 
        std::unordered_map<int, int>& netK);
    bool A_star_routing(Point source, Point sink, int NetId, std::unordered_map<Point,Point,MyHashFunction>& visited_p);
    // return, 0: not find, 1 reach sink, 2 reach tree branch
    int A_star_pin2component_routing(Point source, Point sink, int NetId, std::unordered_map<Point,Point,MyHashFunction>& visited_p,
        std::unordered_map<Point, int, MyHashFunction>& component_map, Point& reach_p);
    int tree2tree_routing(std::priority_queue<std::pair<Point,int>>& p_q, Point b_min, Point b_max, std::unordered_set<int>& source_comp_set,
		std::unordered_set<int>& sink_comp_set, std::vector<std::vector<std::vector<int>>>& comp_grid_map, int NetId, 
		std::unordered_map<Point,Point,MyHashFunction>& visited_p, Point& reach_p);
    int check_segment_profit(Point from, Point to, int NetId);
    TwoPinNet convert_path_to_twopin(Point source, Point sink, std::unordered_map<Point,Point,MyHashFunction>& visited_p);
    TwoPinNet convert_path_to_twopin_t2t(Point& source, Point sink, Point b_min, std::unordered_map<Point,Point,MyHashFunction>& visited_p,
		std::unordered_set<int>& source_comp_set, std::vector<std::vector<std::vector<int>>>& comp_grid_map);
	void set_comp_grid_map(std::vector<std::vector<std::vector<int>>>& comp_grid_map, int netId, Point sink,
        std::unordered_map<Point, int, MyHashFunction>& component_map, Point box_min, Point box_max);
    void set_all_comp_grid_map(std::vector<std::vector<std::vector<int>>>& comp_grid_map, int netId,
        std::unordered_map<Point, int, MyHashFunction>& component_map, Point box_min, Point box_max);
    void add_path_comp_in_pq(std::priority_queue<std::pair<Point,int>>& p_q, std::pair<Point, Point> path, int netId,
        std::unordered_map<Point,Point,MyHashFunction>& visited_p);
    void add_path_comp_in_comp_grid(std::vector<std::vector<std::vector<int>>>& comp_grid_map, 
        std::pair<Point, Point> path, Point box_min, int sink_comp);
    Tree RSMT(std::vector<int> x, std::vector<int> y);
    int row, column, layer;
    int maxCellMove;
    int usedCellMove; //number of moved cells
    SegmentTree* segmentTree;
    int get_edge_cost(Point& from, Point& to);
    std::vector<std::vector<std::unordered_map<int,int>>> cellCount; //for extra demand calculation
    std::vector<MasterCell> masterCells;
    std::vector<Cell> cellInstances;
    std::vector<std::vector<std::vector<Gcell>>> grids;
    std::vector<std::vector<std::unordered_set<int>>> placement;
    std::vector<Net> nets;
    std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, int>>> sameGGrid;
    std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, int>>> adjHGGrid;
    std::unordered_set<int> movedCell; //moved cells
};

#endif
