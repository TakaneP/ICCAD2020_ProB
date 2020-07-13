#include<bits/stdc++.h>
#ifndef PARSER
#define PARSER

class RoutingGraph;

class Parser{
public:
    Parser() = delete;
    Parser(RoutingGraph& main, std::string f): graph(main), file(f)  {};
    void run();
private:
    RoutingGraph& graph;
    std::string file;
};


#endif
