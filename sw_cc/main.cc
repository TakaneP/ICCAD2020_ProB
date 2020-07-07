#include<bits/stdc++.h>
#include"parser.h"
#include"data_structure.h"
#include"segment_tree.h"

using namespace std;

int main(int argc, char* argv[])
{
    RoutingGraph* routingGraph = new RoutingGraph;
    string caseFile = argv[1];
    Parser* parser = new Parser(*routingGraph, caseFile);
    parser->run();
    delete parser;
	
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

    routingGraph->construct_2pin_nets();
    return 0;
}

