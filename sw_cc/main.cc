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
    delete parser;

    //for(int i = 0; i < routingGraph->cellInstances.size(); ++i) {
    //    routingGraph->del_cell(i);
    //}
    for(int i = 0;i < routingGraph->row;i++) {
        for(int j = 0; j < routingGraph->column; j++) {
           for(int k = 0; k < routingGraph->layer;k++) {
                if(k & 1) printf("%d %d %d %d\n", i+1,j+1,k+1,routingGraph->segmentTree->get_remaining_supply(i, i, k, j));
                else printf("%d %d %d %d\n", i+1,j+1,k+1,routingGraph->segmentTree->get_remaining_supply(j, j, k, i));
           }
        }
    }
    return 0;
}

