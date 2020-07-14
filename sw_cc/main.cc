#include<bits/stdc++.h>
#include"parser.h"
#include"data_structure.h"
#include"segment_tree.h"
#include"output.h"

using namespace std;

int main(int argc, char* argv[])
{
    RoutingGraph* routingGraph = new RoutingGraph;
    string caseFile = argv[1];
    Parser* parser = new Parser(*routingGraph, caseFile);
    parser->run();
    delete parser;
    /*
    for(int i = 0; i < routingGraph->cellInstances.size(); ++i) {
        routingGraph->del_cell(i);
    }
    for(int i = 0;i < routingGraph->row;i++) {
        for(int j = 0; j < routingGraph->column; j++) {
            for(int k = 0; k < routingGraph->layer;k++) {
                printf("%d %d %d %d\n", i+1,j+1,k+1,routingGraph->grids[i][j][k].demand);
            }
        }
    }
    */
    /*
	vector<int> x;
	vector<int> y;
	Net& n = routingGraph->nets[0];
	for(auto it =  n.pins.begin(); it != n.pins.end(); ++it) {
	    int cell = it->first;
		x.push_back(routingGraph->cellInstances[cell].x);
		y.push_back(routingGraph->cellInstances[cell].y);
	}
	cout << x.size() << "\n";
	Tree fluteTree = routingGraph->RSMT(x,y);
	cout << fluteTree.length << "\n";
    */
    routingGraph->construct_2pin_nets();
    
    if(argc == 3) {
        string outputFile = argv[2];
        output_file(routingGraph, outputFile);
    }
    return 0;
}

