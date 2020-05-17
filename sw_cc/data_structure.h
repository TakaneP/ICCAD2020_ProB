#include<bits/stdc++.h>
#include<unordered_set>

#ifndef DATASTRUCTURE
#define DATASTRUCTURE

class SegmentTree;

struct Pin{
    Pin() {}
    Pin(int l): layer(l) {}
    int layer;
};

struct MasterCell{
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
    bool operator==(const Point& p2){
        return (this->x==p2.x) && (this->y==p2.y) && (this->z==p2.z);
    }
    int x,y,z;
};

struct Net{
    int minRoutingLayer;
    std::vector<std::pair<int,int>> pins;
    std::vector<std::pair<Point,Point>> routingSegments;
};

struct Cell{
    Cell() {}
    Cell(int mc, bool m, int _x, int _y): masterCell(mc), movable(m), x(_x), y(_y) {}
    int masterCell; 
    bool movable;
    int x,y;
    std::unordered_set<int> connectedNets;    
};

struct Gcell{
    Gcell(): demand(0) {}
    Gcell(int c) : capacity(c), demand(0) {}     
    int get_remaining(void) {return capacity - demand;}
    int capacity;
    int demand;
    std::unordered_set<int> passingNets;
};

class RoutingGraph{
public:
    RoutingGraph();
    ~RoutingGraph();
    void add_cell(int x, int y, int cellIndex);
    void del_cell(int cellIndex);
    void add_cell_demand_into_graph(int x, int y, int MCtype);
    void del_cell_demand_from_graph(int x, int y, int MCtype);
    void add_net_demand_into_graph(int x, int y, int z, int netIndex);
    void del_net_from_graph(int netIndex);
    void del_seg_demand(std::pair<Point,Point> segment, int netIndex);
    void del_seg_demand_from_graph(int x, int y, int z, int netIndex);
    int row, column, layer;
    int maxCellMove;
    SegmentTree* segmentTree;
    std::vector<std::vector<std::unordered_map<int,int>>> cellCount;
    std::vector<MasterCell> masterCells;
    std::vector<Cell> cellInstances;
    std::vector<std::vector<std::vector<Gcell>>> grids;
    std::vector<std::vector<std::unordered_set<int>>> placement;
    std::vector<Net> nets;
    std::vector<std::vector<std::pair<Point,Point>>> routingSegments;
    std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, int>>> sameGGrid;
    std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, int>>> adjHGGrid;
};

#endif
