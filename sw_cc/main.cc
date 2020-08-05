#include<bits/stdc++.h>
#include"parser.h"
#include"data_structure.h"
#include"segment_tree.h"
#include"output.h"

using namespace std;

void print_cell_pin_location(RoutingGraph& routingGraph, int cellIndex, int pinIndex) {
    Cell& cell = routingGraph.cellInstances[cellIndex];
    printf("Cell/Pin: %d/%d , %d %d %d\n", cellIndex, pinIndex, cell.x, cell.y, cell.pins[pinIndex].layer);
}

int main(int argc, char* argv[])
{
    RoutingGraph* routingGraph = new RoutingGraph;
    string caseFile = argv[1];
    Parser* parser = new Parser(*routingGraph, caseFile);
    parser->run();
    delete parser;
    /*vector<int> x;
	vector<int> y;
	Net& n = routingGraph->nets[0];
	for(auto it =  n.pins.begin(); it != n.pins.end(); ++it) {
	    int cell = it->first;
		x.push_back(routingGraph->cellInstances[cell].x);
		y.push_back(routingGraph->cellInstances[cell].y);
	}
    cout << x.size() << "\n";*/
	//Tree fluteTree = routingGraph->RSMT(x,y);
	//cout << fluteTree.length << "\n";
    routingGraph->construct_2pin_nets();
    int prevSize = 0;
    //while(1) {
        routingGraph->move_cells_force();
        //if((int)(routingGraph->movedCell.size() - prevSize) <= 0) break;
        //else prevSize = routingGraph->movedCell.size();
    //}
    if(argc == 3) {
        string outputFile = argv[2];
        output_file(routingGraph, outputFile);
    }
    return 0;
}

