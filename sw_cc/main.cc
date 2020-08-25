#include<bits/stdc++.h>
#include"parser.h"
#include"data_structure.h"
#include"segment_tree.h"
#include"output.h"
#include"global_parameter.h"

using namespace std;

/*********Global Parameter*********/
int failPenalty = 10;
int swapArea = 3;
int moveArea = 5;
int topK = 10;
int layerNormalize = 1;
int mode = 2;
double temperature = 1e2;
/*********Global Parameter*********/


void print_cell_pin_location(RoutingGraph& routingGraph, int cellIndex, int pinIndex) {
    Cell& cell = routingGraph.cellInstances[cellIndex];
    printf("Cell/Pin: %d/%d , %d %d %d\n", cellIndex, pinIndex, cell.x, cell.y, cell.pins[pinIndex].layer);
}

int main(int argc, char* argv[])
{
    time_t start, end, afterPre;
    double diff, onePassTime;
    start = time(NULL);
    random_device rd;
    RoutingGraph* routingGraph = new RoutingGraph(rd);
    string caseFile = argv[1];
    Parser* parser = new Parser(*routingGraph, caseFile);
    parser->run();
    delete parser;
 
    /***********case5********/
    if(routingGraph->row == 104 && routingGraph->column == 103 && routingGraph->layer == 16) {
        layerNormalize = routingGraph->layer;
    }
    /***********case5********/

    routingGraph->construct_2pin_nets();
    int prevSize = 0;
    int count = 0;
    // 2 pin reroute
    int wire_length=0;
    mode = 2;
    routingGraph->wirelength_driven_move(wire_length);
    end = time(NULL);
    afterPre = end;
    diff = difftime(end, start);
    bool firstTime = true;
    if(diff < 3000.0) {
        while(1) {
            //routingGraph->move_cells_force();
            int wire_length=0;
            mode = count % 2;
            routingGraph->wirelength_driven_move(wire_length);
            //cout << "####size: " << routingGraph->movedCell.size() << endl;
            //cout << "####imporve " << wire_length << endl;
            if(wire_length <= 0) break;
            // if((int)(routingGraph->movedCell.size() - prevSize) <= 0) break;
            // else prevSize = routingGraph->movedCell.size();
            count++;
            temperature = (temperature>10)? temperature/10:1;
            for(auto& net : routingGraph->nets) {
                net.remove_dangling_wire(routingGraph->grids);
            }

            if(firstTime) {
                firstTime = false;
                end = time(NULL);
                onePassTime = difftime(end, afterPre);
            }
            end = time(NULL);
            diff = difftime(end, start);
            if((3600.0 - diff) <= onePassTime) break;
        }
    }
    if(argc == 3) {
        string outputFile = argv[2];
        output_file(routingGraph, outputFile);
    }
    return 0;
}

