#include<bits/stdc++.h>
#include "./parser.h"
#include "data_structure.h"
#include "segment_tree.h"
using namespace std;

int get_postfix_int(string& s, int start) {
    return stoi(s.substr(start)) -1;
}

void Parser::run(void) {
    ifstream fin;
	fin.open(file);
    string type;
    int value;
    //Cell max move
    fin >> type >> graph.maxCellMove;
    //Grid boundaries
    fin >> type >> value >> value >> graph.row;
    fin >> graph.column;
    double segmentTreeRowNodeNum_D = log2(graph.row);
    double segmentTreeColNodeNum_D = log2(graph.column);
    segmentTreeRowNodeNum_D = pow(2,ceil(segmentTreeRowNodeNum_D)+1);
    segmentTreeColNodeNum_D = pow(2,ceil(segmentTreeColNodeNum_D)+1);
    int segmentTreeRowNodeNum = (int) segmentTreeRowNodeNum_D;
    int segmentTreeColNodeNum = (int) segmentTreeColNodeNum_D;
    //Layer
    fin >> type >> graph.layer;
    vector<int> layerLimit(graph.layer);
    graph.segmentTree->node.resize(graph.layer);
    for(int i = 0; i < graph.layer; i++) {
        fin >> type >> type >> type >> type >> layerLimit[i];
        if(i & 1) {
            graph.segmentTree->node[i].resize(graph.column);
            for(int j = 0; j < graph.column; j++) {
                graph.segmentTree->node[i][j].resize(segmentTreeRowNodeNum);
            }
        }
	    else {
            graph.segmentTree->node[i].resize(graph.row);
            for(int j = 0; j < graph.row; j++) {
                graph.segmentTree->node[i][j].resize(segmentTreeColNodeNum);
            }
        }
    }
    //Initial each layer's routing supply
    for(int i = 0; i < graph.row; i++) {
        graph.grids.push_back(vector<vector<Gcell>>());
        graph.cellCount.push_back(vector<unordered_map<int,int>>());
        graph.placement.push_back(vector<unordered_set<int>>());
        for(int j = 0; j < graph.column; j++) {
            graph.grids[i].push_back(vector<Gcell>());
            graph.cellCount[i].push_back(unordered_map<int,int>());
            graph.placement[i].push_back(unordered_set<int>());
            for(int k = 0; k < graph.layer; k++) {
                graph.grids[i][j].push_back(Gcell(layerLimit[k]));
	        }
        }
    }
    //Non default supply
    fin >> type >> value;
    for(int i = 0; i < value; ++i) {
        int x,y,z,offset;
	    fin >> x >> y >> z >> offset;
	    graph.grids[x-1][y-1][z-1].capacity += offset;
    }
    //Master Cell
    int pinNum, blkNum;
    fin >> type >> value;
    graph.masterCells.resize(value);
    for(int i = 0; i < value; ++i) {
        fin >> type >> type >> pinNum >> blkNum;
        graph.masterCells[i].pins.resize(pinNum);
        for(int j = 0;j < pinNum; j++) {
            fin >> type >> type >> type;
            int pinLayer = get_postfix_int(type, 1);
            graph.masterCells[i].pins[j] = Pin(pinLayer);
        }
        for(int j = 0;j < blkNum; j++) {
            int demand, blkLayer;
            fin >> type >> type >> type >> demand;
            blkLayer = get_postfix_int(type, 1);
            graph.masterCells[i].blockages[blkLayer] += demand;
        }
    }
    //Extra Demand
    fin >> type >> value;
    for(int i = 0;i < value; ++i) {
        bool op;
        int MC1, MC2, layerIndex, extraDemand;
        fin >> type;
        op = (type[0] == 'a'); //same : 0, adj: 1
        fin >> type;
        MC1 = get_postfix_int(type, 2);
        fin >> type;
        MC2 = get_postfix_int(type, 2);
        fin >> type;
        layerIndex = get_postfix_int(type, 1);
        fin >> extraDemand;
        if(!op) {
            graph.sameGGrid[MC1][MC2][layerIndex] = extraDemand;
            graph.sameGGrid[MC2][MC1][layerIndex] = extraDemand;
        }
        else {
            graph.adjHGGrid[MC1][MC2][layerIndex] = extraDemand;
            graph.adjHGGrid[MC2][MC1][layerIndex] = extraDemand;
        }
    }
    //Cell Instances
    fin >> type >> value;
    graph.cellInstances.resize(value);
    for(int i = 0;i < value; ++i) {
        int rowPos, colPos, MCtype;
        fin >> type >> type >> type;
        MCtype = get_postfix_int(type,2); //MC
        fin >> rowPos >> colPos >> type;
        bool mov = (type[0] == 'M');
        MasterCell& m = graph.masterCells[MCtype]; 
        graph.cellInstances[i] = Cell(m.pins, m.blockages, mov, --rowPos, --colPos, MCtype);
        graph.add_cell(rowPos, colPos, i);
    }
    //Net
    fin >> type >> value;
    graph.nets.resize(value);
    for(int i = 0; i < value; ++i) {
        graph.nets[i].netId = i;
        int pinNum;
        fin >> type >> type >> pinNum >> type;
        graph.nets[i].minRoutingLayer = (type[0] == 'N')? 0:get_postfix_int(type,1);
        for(int j = 0; j < pinNum; j++) {
            fin >> type >> type;
            size_t found = type.find("/");
            string cellIns = type.substr(0,found);
            string pinName = type.substr(found+1);
            int cellIndex = get_postfix_int(cellIns, 1);
            int pinIndex = get_postfix_int(pinName, 1);
            graph.nets[i].pins.push_back({cellIndex, pinIndex});
            Cell& cell = graph.cellInstances[cellIndex];
            cell.pins[pinIndex].connectedNet = i;
            graph.nets[i].add_net_demand_into_graph(cell.x, cell.y, cell.pins[pinIndex].layer, graph.grids);
        }
    }
    //Routing segments
    fin >> type >> value;
    for(int i = 0;i < value; ++i) {
        int startRow, startColumn, startLayer, endRow, endColumn, endLayer;
        fin >> startRow >> startColumn >> startLayer >> endRow >> endColumn >> endLayer >> type;
        int netIndex = get_postfix_int(type,1);
        startRow--;
        startColumn--;
        startLayer--;
        endRow--;
        endColumn--;
        endLayer--;
        if(startRow > endRow) swap(startRow, endRow);
        else if(startColumn > endColumn) swap(startColumn, endColumn);
        else if(startLayer > endLayer) swap(startLayer, endLayer);
        graph.nets[netIndex].routingSegments.push_back({Point( startRow, startColumn, startLayer), 
                                                        Point(endRow, endColumn, endLayer)});
        if(startRow != endRow) {
            for(int j = startRow; j <= endRow; j++) {
                graph.nets[netIndex].add_net_demand_into_graph(j, startColumn, startLayer, graph.grids);
            }
        }
        else if(startColumn != endColumn) {
            for(int j = startColumn; j <= endColumn; j++) {
                graph.nets[netIndex].add_net_demand_into_graph(startRow, j, startLayer, graph.grids);
            }
        }
        else if(startLayer != endLayer) {
            for(int j = startLayer; j <= endLayer; j++) {
                graph.nets[netIndex].add_net_demand_into_graph(startRow, startColumn, j, graph.grids);
            }
        }
    }
    graph.segmentTree->build_ini();
	fin.close();
}
