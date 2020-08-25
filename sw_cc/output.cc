#include<bits/stdc++.h>
#include"output.h"
#include"data_structure.h"

using namespace std;

void DFS(vector<pair<Point,Point>>& netSegments, map<tuple<int,int,int>, bool>& visit, unordered_map<Point,TreeNode,MyHashFunction>& branchNodes,const Point& cur, int& numRoutes) {
    TreeNode& treeNode = branchNodes[cur];
    visit[{cur.x, cur.y, cur.z}] = true; 
    for(auto it = treeNode.neighbors.begin(); it != treeNode.neighbors.end(); ++it) {
        const Point& nxt = it->first;
        TwoPinNet& twopinNet = it->second;
        if(visit[{nxt.x, nxt.y, nxt.z}] == false) {
            for(auto _it = twopinNet.paths.begin(); _it != twopinNet.paths.end(); ++_it) {
                netSegments.push_back(*_it);
                numRoutes++;
            }
            DFS(netSegments, visit, branchNodes, nxt, numRoutes);
        }
    }
}

void output_file(RoutingGraph* routingGraph, string outputFile) {
    ofstream fout;
    fout.open(outputFile);
    fout << "NumMovedCellInst " << routingGraph->movedCell.size() << "\n";
    for(auto& cell : routingGraph->movedCell) 
        fout << "CellInst " << "C" << cell + 1 << " " << routingGraph->cellInstances[cell].x + 1 << " " << routingGraph->cellInstances[cell].y + 1 << "\n";
    int numRoutes = 0;
    vector<vector<pair<Point,Point>>> storeSegments(routingGraph->nets.size());
    for(int i = 0; i < routingGraph->nets.size(); ++i) {
        Net& net = routingGraph->nets[i];
        map<tuple<int,int,int>, bool> visit;
        if(!net.branch_nodes.empty()) {
            const Point& p = net.branch_nodes.begin()->first;
            DFS(storeSegments[i], visit, net.branch_nodes, p, numRoutes);
        }
        if(net.minRoutingLayer > 0) {
            map<pair<int, int>, set<int>> insertVias;
            for(auto it = net.pins.begin(); it != net.pins.end(); ++it) {
                int cellIdx = it->first;
                int pinIdx = it->second;
                Cell& cell = routingGraph->cellInstances[cellIdx];
                Pin& pin = cell.pins[pinIdx];
                if(pin.pseudo && pin.actualPinLayer < net.minRoutingLayer)
                    insertVias[{cell.x, cell.y}].insert(pin.actualPinLayer);
            }
            for(auto it = insertVias.begin(); it != insertVias.end(); ++it) {
                if(it->second.empty()) continue;
                int minLayer = *(it->second.begin());
                int x = it->first.first;
                int y = it->first.second;
                storeSegments[i].push_back({Point(x, y, minLayer), Point(x, y, net.minRoutingLayer)});
                numRoutes++;
            }
        }
    }
    fout << "NumRoutes " << numRoutes << "\n";
    for(int i = 0; i < storeSegments.size(); ++i) {
        vector<pair<Point, Point>>& v = storeSegments[i];
        for(auto it = v.begin(); it != v.end(); ++it) {
            Point p1 = it->first;
            Point p2 = it->second;
            fout << p1.x+1 << " " << p1.y+1 << " " << p1.z+1 << " " << p2.x+1 << " " << p2.y+1 << " " << p2.z+1 << " N" << i+1 << "\n";
        }
    }
}

int get_output_wirelength(RoutingGraph* routingGraph, vector<vector<pair<Point,Point>>>& storeSegments, int& routesNumber) {
    int numRoutes = 0;
    //vector<vector<pair<Point,Point>>> storeSegments(routingGraph->nets.size());
    for(int i = 0; i < routingGraph->nets.size(); ++i) {
        Net& net = routingGraph->nets[i];
        map<tuple<int,int,int>, bool> visit;
        if(!net.branch_nodes.empty()) {
            const Point& p = net.branch_nodes.begin()->first;
            DFS(storeSegments[i], visit, net.branch_nodes, p, numRoutes);
        }
        if(net.minRoutingLayer > 0) {
            map<pair<int, int>, set<int>> insertVias;
            for(auto it = net.pins.begin(); it != net.pins.end(); ++it) {
                int cellIdx = it->first;
                int pinIdx = it->second;
                Cell& cell = routingGraph->cellInstances[cellIdx];
                Pin& pin = cell.pins[pinIdx];
                if(pin.pseudo && pin.actualPinLayer < net.minRoutingLayer)
                    insertVias[{cell.x, cell.y}].insert(pin.actualPinLayer);
            }
            for(auto it = insertVias.begin(); it != insertVias.end(); ++it) {
                if(it->second.empty()) continue;
                int minLayer = *(it->second.begin());
                int x = it->first.first;
                int y = it->first.second;
                storeSegments[i].push_back({Point(x, y, minLayer), Point(x, y, net.minRoutingLayer)});
                numRoutes++;
            }
        }
    }
    int wirelength = 0;
    for(int i = 0; i < storeSegments.size(); ++i) {
        vector<pair<Point, Point>>& v = storeSegments[i];
        for(auto it = v.begin(); it != v.end(); ++it) {
            Point p1 = it->first;
            Point p2 = it->second;
            wirelength += ADIFF(p1.x, p2.x) + ADIFF(p1.y, p2.y) + ADIFF(p1.z, p2.z);
        }
    }
    routesNumber = numRoutes;
    return wirelength;
}

void output_file_with_segments(RoutingGraph* routingGraph, string outputFile, vector<vector<pair<Point,Point>>>& storeSegments, int& routesNumber) {
    ofstream fout;
    fout.open(outputFile);
    fout << "NumMovedCellInst " << routingGraph->movedCell.size() << "\n";
    for(auto& cell : routingGraph->movedCell)
        fout << "CellInst " << "C" << cell + 1 << " " << routingGraph->cellInstances[cell].x + 1 << " " << routingGraph->cellInstances[cell].y + 1 << "\n";
    fout << "NumRoutes " << routesNumber << "\n";
    for(int i = 0; i < storeSegments.size(); ++i) {
        vector<pair<Point, Point>>& v = storeSegments[i];
        for(auto it = v.begin(); it != v.end(); ++it) {
            Point p1 = it->first;
            Point p2 = it->second;
            fout << p1.x+1 << " " << p1.y+1 << " " << p1.z+1 << " " << p2.x+1 << " " << p2.y+1 << " " << p2.z+1 << " N" << i+1 << "\n";
        }
    }
}
