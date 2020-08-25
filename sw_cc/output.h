#ifndef _OUTPUTFILE_
#define _OUTPUTFILE_
#include<bits/stdc++.h>
#include"data_structure.h"

class RoutingGraph;

void output_file(RoutingGraph* routingGraph, std::string outputFile);
int get_output_wirelength(RoutingGraph* routingGraph, std::vector<std::vector<std::pair<Point,Point>>>& storeSegments, int& routesNumber);
void output_file_with_segments(RoutingGraph* routingGraph, std::string outputFile, std::vector<std::vector<std::pair<Point,Point>>>& storeSegments, int& routesNumber);
#endif
