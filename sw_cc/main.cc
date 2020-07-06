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

    routingGraph->construct_2pin_nets();
    return 0;
}

