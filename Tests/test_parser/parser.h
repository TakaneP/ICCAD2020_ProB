#include<unordered_set>
#include<unordered_map>
#include<string>
#include<vector>
#include<utility>
#include<tuple>

#ifndef PARSER
#define PARSER

struct Blockage{
    Blockage() {}
    Blockage(int l, int d): layer(l), demand(d) {}
    int layer;
    int demand;
};

struct Pin{
    Pin() {}
    Pin(int l): layer(l) {}
    int layer;
};

struct MasterCell{
    std::unordered_map<std::string, Pin> pins;
    std::unordered_map<std::string, Blockage> blockages;
};

struct Net{
    int minRoutingLayer;
    std::vector<std::pair<std::string,std::string>> pins;
};

struct Cell{
    Cell() {}
    Cell(std::string mc, bool m): masterCellName(mc), movable(m) {}
    
    std::string masterCellName; 
    bool movable;
    int x,y;    
};

class Gcell{
public:
    Gcell(): demand(0) {}
    Gcell(int c) : capacity(c), demand(0) {} 
    
    int get_capacity(void) {return capacity;}
    int get_demand(void) {return demand;}
    std::unordered_set<std::string>& get_passingNets(void) {return passingNets;}
    
    void add_net_demand(void) {demand += passingNets.size();}
    void set_capacity(int c) {capacity = c;}
    void adjust_capacity(int offset) {capacity += offset;}
    void set_demand(int d) {demand = d;}
    void adjust_demand(int offset) {demand += offset;}
    void insert_net(std::string n) {passingNets.insert(n);} 
private:
    int capacity;
    int demand;
    std::unordered_set<std::string> passingNets;
};

struct Point{
    Point() {}
    Point(int i, int j, int k): x(i), y(j), z(k) {}
    int x,y,z;
};

class RoutingGraph{
public: 
    RoutingGraph(int r, int c, int l, int m);
    int get_row(void) {return row;}
    int get_column(void) {return column;}
    int get_layer(void) {return layer;}
    int get_maxCellMove(void) {return maxCellMove;}

    MasterCell& get_masterCell(std::string s) {return masterCells[s];}
    int& get_layer(std::string s) {return layers[s];}
    Cell& get_cellInstance(std::string s) {return cellInstances[s];}
    Net& get_net(std::string s) {return nets[s];}
    std::vector<std::pair<Point,Point>>& get_routingSegments(std::string s) {return routingSegments[s];}

    Gcell& get_Gcell(int i, int j, int k) {return grids[i][j][k];}
    
    std::vector<std::vector<std::vector<Gcell>>>& get_grids(void) {return grids;}
    std::unordered_map<std::string, int>& get_layers(void) {return layers;}
    std::vector<std::vector<std::unordered_set<std::string>>>& get_placement(void) {return placement;}
    std::unordered_map<std::string, MasterCell>& get_masterCells(void) {return masterCells;}
    std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<int, int>>>& get_sameGGrid(void) {return sameGGrid;}
    std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<int, int>>>& get_adjHGGrid(void) {return adjHGGrid;}
    std::unordered_map<std::string, Cell>& get_cellInstances(void) {return cellInstances;}
    std::unordered_map<std::string, Net>& get_nets(void) {return nets;}
    std::unordered_map<std::string, std::vector<std::pair<Point,Point>>>& get_routingSegments(void) {return routingSegments;}
private:
    int row, column, layer;
    int maxCellMove;
    std::unordered_map<std::string, int> layers;
    std::unordered_map<std::string, MasterCell> masterCells;
    std::unordered_map<std::string, Cell> cellInstances;
    std::vector<std::vector<std::vector<Gcell>>> grids;
    std::vector<std::vector<std::unordered_set<std::string>>> placement;
    std::unordered_map<std::string, Net> nets;
    std::unordered_map<std::string, std::vector<std::pair<Point,Point>>> routingSegments;
    std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<int, int>>> sameGGrid;
    std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<int, int>>> adjHGGrid;
};

#endif
