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
	
	vector<int> x;
	vector<int> y;
	Net& n = routingGraph->nets[0];
	for(auto it =  n.pins.begin(); it != n.pins.end(); ++it) {
	    int cell = it->first;
		x.push_back(routingGraph->cellInstances[cell].x);
		y.push_back(routingGraph->cellInstances[cell].y);
	}
	cout << x.size() << "\n";
	//Tree fluteTree = routingGraph->RSMT(x,y);
	//cout << fluteTree.length << "\n";

    routingGraph->construct_2pin_nets();
    routingGraph->move_cells_force();
    /*Point insert_p(1,2,1);
    routingGraph->del_cell_neighbor(7);
    routingGraph->add_cell(2,2,7);
    unordered_map<Point, int, MyHashFunction> component_map;
    routingGraph->nets[1].set_point_component(component_map);
    unordered_map<Point,Point,MyHashFunction> visited_p;
    Point source(2,2,0), sink1(4,1,0), sink2(1,3,0), reach_p;
    int find_flg = routingGraph->A_star_pin2component_routing(source, sink1, 1, visited_p, component_map, reach_p);
    cout << "find flg: " << find_flg << " reach_p: " << reach_p << endl;   
    TwoPinNet two_pin = routingGraph->convert_path_to_twopin(source, reach_p, visited_p);

    if(routingGraph->nets[1].branch_nodes.find(reach_p) == routingGraph->nets[1].branch_nodes.end())
        routingGraph->nets[1].insert_steiner_point(reach_p);
    if(source == reach_p)
        routingGraph->nets[1].branch_nodes[reach_p].node.type = 1;
    else {
        routingGraph->nets[1].branch_nodes[source].neighbors.emplace_back(reach_p, two_pin);
        routingGraph->nets[1].branch_nodes[reach_p].neighbors.emplace_back(source, two_pin_reverse(two_pin));
    }
    
    cout << "after find\n";
    for(auto& node : routingGraph->nets[1].branch_nodes) {
        cout << node.first << endl;
        for(auto& neighbor : node.second.neighbors) {
            cout << "\t" << neighbor.first << endl;
            for(auto& path : neighbor.second.paths) {
                cout << "\t\t" << path.first << "->" << path.second << endl;
            }
        }
    }
    unordered_map<int, int> netK;
    netK[1]++;
    routingGraph->del_cell_last_k_neighbor(7, netK);
    routingGraph->nets[1].clear_steiner_point(reach_p, routingGraph->grids);
    cout << "after merge\n";
    for(auto& node : routingGraph->nets[1].branch_nodes) {
        cout << node.first << endl;
        for(auto& neighbor : node.second.neighbors) {
            cout << "\t" << neighbor.first << endl;
            for(auto& path : neighbor.second.paths) {
                cout << "\t\t" << path.first << "->" << path.second << endl;
            }
        }
    }*/
    if(argc == 3) {
        string outputFile = argv[2];
        output_file(routingGraph, outputFile);
    }
    return 0;
}

