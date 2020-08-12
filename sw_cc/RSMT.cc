#include<bits/stdc++.h>
#include"data_structure.h"
#include"RSMT.h"

using namespace std;

RSMT::RSMT(RoutingGraph& g): routingGraph(g) {};

void RSMT::RSMT_reroute() {
    //for all nets
    for(auto& net : routingGraph.nets) {
        //DFS release routing resource, record original length    
        if(net.branch_nodes.empty()) continue;
        int length = 0;
        unordered_map<Point, bool, MyHashFunction> visited;
        stack<Point> candidate;
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
                net.del_twoPinNet_from_graph(twoPinNet, routingGraph.grids);
                twoPinNet.update_wire_length();
                length += twoPinNet.wire_length;
                if(visited[neighbor] != true) {
                    visited[neighbor] = true;
                    candidate.push(neighbor);
                }
            }
        }
        //get flute tree
        vector<int> x;
	    vector<int> y;
	    for(auto it =  net.pins.begin(); it != net.pins.end(); ++it) {
	        int cell = it->first;
	    	x.push_back(routingGraph.cellInstances[cell].x);
	    	y.push_back(routingGraph.cellInstances[cell].y);
	    }
	    Tree fluteTree = routingGraph.RSMT(x,y);
        for(int i = 0; i < 2*fluteTree.deg-2; ++i) {
            printf("%d %d\n", fluteTree.branch[i].x+1, fluteTree.branch[i].y+1);
            printf("%d %d\n\n", fluteTree.branch[fluteTree.branch[i].n].x+1,
               fluteTree.branch[fluteTree.branch[i].n].y+1);
        }
        //get pin position
        unordered_map<int, unordered_map<int, int>> minLayer;
        unordered_map<int, unordered_map<int, int>> maxLayer;
        for(auto it = net.branch_nodes.begin(); it != net.branch_nodes.end(); ++it) {
            const Point& p = it->first;
            TreeNode& treeNode = it->second;
            if(treeNode.node.type == 2 || treeNode.node.type == 1) {
                minLayer[p.x][p.y] = min(minLayer[p.x][p.y], p.z);
                maxLayer[p.x][p.y] = max(maxLayer[p.x][p.y], p.z);
            }
        }
        //Fake temp routing tree
            

    }
 
}
