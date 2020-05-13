#include<bits/stdc++.h>
#include"parser.h"
#include"data_structure.h"

using namespace std;

int main(int argc, char* argv[])
{
    RoutingGraph* routingGraph = new RoutingGraph;
    string caseFile = argv[1];
    Parser* parser = new Parser(*routingGraph, caseFile);
    parser->run();
    //for testing
    /*for(int i = 0;i < routingGraph->row;i++) {
        for(int j = 0; j < routingGraph->column; j++) {
	    for(int k = 0; k < routingGraph->layer;k++) {
                printf("%d %d %d %d\n", i+1,j+1,k+1,routingGraph->grids[i][j][k].demand);
	    }
	}
    }*/
    return 0;
}

