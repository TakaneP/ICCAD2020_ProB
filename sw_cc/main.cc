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
int predictRegionSize = 1;
bool case1 = 0;
bool case2 = 0;
bool case3 = 0;
bool case4 = 0;
bool case5 = 0;
bool case6 = 0;
bool openSA = 0;
/*********Global Parameter*********/


void print_cell_pin_location(RoutingGraph& routingGraph, int cellIndex, int pinIndex) {
    Cell& cell = routingGraph.cellInstances[cellIndex];
    printf("Cell/Pin: %d/%d , %d %d %d\n", cellIndex, pinIndex, cell.x, cell.y, cell.pins[pinIndex].layer);
}

int main(int argc, char* argv[])
{
    time_t start, end, afterPre, roundStart;
    double diff, onePassTime, firstRoundTime;
    start = time(NULL);
    int minWireLength = INT32_MAX;
    bool firstRound = true;
    while(1) {
        temperature = 1e2;
        roundStart = time(NULL);
        random_device rd;
        RoutingGraph* routingGraph = new RoutingGraph(rd);
        string caseFile = argv[1];
        Parser* parser = new Parser(*routingGraph, caseFile);
        parser->run();
        delete parser;
 
        /***********case1*******/
        if(routingGraph->row == 5 && routingGraph->column == 5 && routingGraph->layer == 3) {
            case1 = true;
            predictRegionSize = 100;
        }
        /***********case1********/

        /***********case2*******/
        if(routingGraph->row == 4 && routingGraph->column == 4 && routingGraph->layer == 3) {
            case2 = true;
        }
        /***********case2********/

        /***********case3*******/
        if(routingGraph->row == 27 && routingGraph->column == 33 && routingGraph->layer == 7) {
            case3 = true;
        }
        /***********case3********/

        /***********case4*******/
        if(routingGraph->row == 277 && routingGraph->column == 277 && routingGraph->layer == 12) {
            case4 = true;
        }
        /***********case4********/

        /***********case5********/
        if(routingGraph->row == 104 && routingGraph->column == 103 && routingGraph->layer == 16) {
            case5 = true;
            layerNormalize = routingGraph->layer;
        }
        /***********case5********/

        /***********case6*******/
        if(routingGraph->row == 237 && routingGraph->column == 236 && routingGraph->layer == 16) {
            case6 = true;
        }
        /***********case6********/

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
                if((3600.0 - diff) <= onePassTime*2.0 || (3600.0 - diff) <= 600) break;
            }
        }
        int routesNum = 0;
        vector<vector<pair<Point, Point>>> segments(routingGraph->nets.size());
        int currentWL = get_output_wirelength(routingGraph, segments, routesNum);
        if(currentWL < minWireLength) {
            if(argc == 3) {
                string outputFile = argv[2];
                output_file_with_segments(routingGraph, outputFile, segments, routesNum);
            }
            minWireLength = currentWL;
        }
        if(firstRound) {
            firstRound = false;
            end = time(NULL);
            firstRoundTime = difftime(end, roundStart);
        }
        end = time(NULL);
        diff = difftime(end, start);
        //cout << currentWL << " " << 3600-diff << " " << firstRoundTime << endl;
        if((3600 - diff) <= firstRoundTime*2.0) break;
        if(case1 || case2) return 0;
        delete routingGraph;
        if(openSA == false) openSA = true;
    }
    return 0;
}

