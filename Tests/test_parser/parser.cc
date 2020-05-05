#include<iostream>
#include<fstream>
#include "./parser.h"
using namespace std;

RoutingGraph::RoutingGraph(int r, int c, int l, int m): row(r), column(c), layer(l), maxCellMove(m) {
    grids.resize(r);
    placement.resize(r);
    for(int i = 0;i < r;i++) {
        grids[i].resize(c);
	placement[i].resize(c);
	for(int j = 0; j < c;j++) {
	    grids[i][j].resize(l);
	}
    }
}

RoutingGraph* parser(RoutingGraph* routingGraph, string file)
{
    ifstream fin(file);
    string type;
    int maxCellMove;
    int rowBegin, rowEnd, columnBegin, columnEnd;
    int numLayer;
    fin >> type >> maxCellMove;
    fin >> type >> rowBegin >> columnBegin >> rowEnd >> columnEnd;
    fin >> type >> numLayer;
    int rowNum = rowEnd - rowBegin + 1;
    int columnNum = columnEnd - columnBegin + 1;
    routingGraph = new RoutingGraph(rowNum, columnNum, numLayer, maxCellMove);
    vector<vector<vector<Gcell>>>& grids = routingGraph->get_grids();
    
    unordered_map<string,int>& layers = routingGraph->get_layers();
    int address;
    for(int i = 0; i < numLayer ; ++i) {
    	string layerName, direction;
	int index, supply;
	fin >> type >> layerName >> index >> direction >> supply;
	layers[layerName] = --index;
        
	for(int j = 0; j < rowNum; ++j) {
	    for(int k = 0; k < columnNum; ++k) {
	        grids[j][k][index].set_capacity(supply);
	    }
	}
    }
    int numNonDefaultSupplyGGrid;
    fin >> type >> numNonDefaultSupplyGGrid;
    for(int i = 0; i < numNonDefaultSupplyGGrid; ++i) {
        int x,y,z,offset;
	fin >> x >> y >> z >> offset;
	grids[x-1][y-1][z-1].adjust_capacity(offset);
    }
    int numMasterCell;
    string masterCellName;
    int pinNum, blkNum;
    fin >> type >> numMasterCell;
    for(int i = 0; i < numMasterCell; ++i) {
        fin >> type >> masterCellName >> pinNum >> blkNum;
	MasterCell& masterCell = routingGraph->get_masterCell(masterCellName);
	unordered_map<string, Pin> pins;
        unordered_map<string, Blockage> blks;	
	for(int j = 0;j < pinNum; j++) {
	    string pinName, pinLayer;
	    fin >> type >> pinName >> pinLayer;
	    pins[pinName] = Pin(routingGraph->get_layer(pinLayer));
	}
	for(int j = 0;j < blkNum; j++) {
	    string blkName, blkLayer;
	    int demand;
	    fin >> type >> blkName >> blkLayer >> demand;
	    blks[blkName] = Blockage(routingGraph->get_layer(blkLayer), demand);
	}
	masterCell.pins = pins;
	masterCell.blockages = blks;
    }
   
    int extraDemandNum;
    fin >> type >> extraDemandNum;
    unordered_map<string, unordered_map<string, unordered_map<int, int>>>& sameGGrid = routingGraph->get_sameGGrid();
    unordered_map<string, unordered_map<string, unordered_map<int, int>>>& adjHGGrid = routingGraph->get_adjHGGrid();
    for(int i = 0;i < extraDemandNum; ++i) {
        string MC1, MC2, layerName;
	int extraDemand;
	fin >> type >> MC1 >> MC2 >> layerName >> extraDemand;
	int layerIndex = routingGraph->get_layer(layerName);
	if(type == "sameGGrid") sameGGrid[MC1][MC2][layerIndex] = extraDemand;
	else if(type == "adjHGGrid") {
	    adjHGGrid[MC1][MC2][layerIndex] = extraDemand;
	}
    }
    int numCellIns;
    fin >> type >> numCellIns;
    vector<vector<unordered_set<string>>>& placement = routingGraph->get_placement();
    for(int i = 0;i < numCellIns; ++i) {
        string cellName, masterCellName, movable;
	int rowPos, colPos;
	fin >> type >> cellName >> masterCellName >> rowPos >> colPos >> movable;
	placement[rowPos-1][colPos-1].insert(cellName);
	Cell& cell = routingGraph->get_cellInstance(cellName);
	cell.masterCellName = masterCellName;
	cell.movable = (movable == "Movable")? true:false;
	cell.x = rowPos-1;
	cell.y = colPos-1;
    }
    int numNets;
    fin >> type >> numNets;
    for(int i = 0; i < numNets; ++i) {
        string netName, minLayerCons;
	int pinNum;
	fin >> type >> netName >> pinNum >> minLayerCons;
	Net& net = routingGraph->get_net(netName);
	net.minRoutingLayer = (minLayerCons == "NoCstr")? -1:routingGraph->get_layer(minLayerCons);
	for(int j = 0; j < pinNum; j++) {
	    string component;
	    fin >> type >> component;
	    size_t found = component.find("/");
	    string cellIns = component.substr(0,found);
	    string pinName = component.substr(found+1);
	    net.pins.push_back(make_pair(cellIns, pinName));
	}
    }
    int numRoutes;
    fin >> type >> numRoutes;
    for(int i = 0;i < numRoutes; ++i) {
        int startRow, startColumn, startLayer, endRow, endColumn, endLayer;
	string netName;
	fin >> startRow >> startColumn >> startLayer >> endRow >> endColumn >> endLayer >> netName;
	vector<pair<Point,Point>>& netSegments = routingGraph->get_routingSegments(netName);
	netSegments.push_back(make_pair(Point(startRow-1, startColumn-1, startLayer-1), Point(endRow-1, endColumn-1, endLayer-1)));
    }
    fin.close();
   
    return routingGraph;
}

void debug_print(RoutingGraph* routingGraph)
{
    int row = routingGraph->get_row();
    int column = routingGraph->get_column();
    int layer = routingGraph->get_layer();
    printf("%d\n", routingGraph->get_maxCellMove());
    printf("%d %d %d\n", row, column, layer);
    std::unordered_map<std::string, int>& layers = routingGraph->get_layers();
    for(auto it = layers.begin(); it != layers.end(); ++it) {
	printf("%s %d\n", it->first.c_str(), it->second+1);
    }
    std::vector<std::vector<std::vector<Gcell>>>& grids = routingGraph->get_grids();
    for(int i = 0; i < row; i++) {
        for(int j = 0; j < column; j++) {
	    for(int k = 0; k < layer; k++) {
	        printf("%d %d %d supply %d\n", i+1,j+1,k+1,grids[i][j][k].get_capacity());
	    }
	}
    }
    unordered_map<string, MasterCell> masterCells = routingGraph->get_masterCells();
    for(auto it = masterCells.begin(); it != masterCells.end(); ++it) {
        MasterCell& cur = it->second;
	for(auto iit = cur.pins.begin(); iit != cur.pins.end(); ++iit) {
	    printf("%s %s %d\n", it->first.c_str(), iit->first.c_str(), iit->second.layer+1);
	}
	for(auto iit = cur.blockages.begin(); iit != cur.blockages.end(); ++iit) {
	    printf("%s %s %d %d\n", it->first.c_str(), iit->first.c_str(), iit->second.layer+1,iit->second.demand);
	}
    }
    std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<int, int>>>& sameGGrid = routingGraph->get_sameGGrid();
    std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<int, int>>>& adjHGGrid = routingGraph->get_adjHGGrid();
    for(auto it = sameGGrid.begin(); it != sameGGrid.end(); ++it) {
        for(auto iit = it->second.begin(); iit != it->second.end(); ++iit) {
	    for(auto _it = iit->second.begin(); _it != iit->second.end(); ++_it) {
	        printf("sameGGrid %s %s %d %d\n", it->first.c_str(), iit->first.c_str(), _it->first+1, _it->second);
	    }
	}
    }
    for(auto it = adjHGGrid.begin(); it != adjHGGrid.end(); ++it) {
        for(auto iit = it->second.begin(); iit != it->second.end(); ++iit) {
            for(auto _it = iit->second.begin(); _it != iit->second.end(); ++_it) {
                printf("adjHGGrid %s %s %d %d\n", it->first.c_str(), iit->first.c_str(), _it->first+1, _it->second);
            }
        }
    }
    std::unordered_map<std::string, Cell>& cellInstances = routingGraph->get_cellInstances();
    std::vector<std::vector<std::unordered_set<std::string>>>& placement = routingGraph->get_placement();
    for(int i = 0; i < row; ++i) {
        for(int j = 0;j < column; ++j) {
	    if(!placement[i][j].empty()) {
	        for(auto it = placement[i][j].begin(); it != placement[i][j].end(); ++it) {
		    Cell& cell = cellInstances[*it];
		    printf("%d %d %s %s %s\n",i+1,j+1, it->c_str(),cell.masterCellName.c_str(), (cell.movable == true)? "Movable":"Fixed");
		}
	    }
	}
    }
    std::unordered_map<std::string, Net>& nets = routingGraph->get_nets();
    for(auto it = nets.begin(); it != nets.end(); ++it) {
        Net& net = it->second;
	for(auto _it = net.pins.begin(); _it != net.pins.end(); ++_it) {
	    printf("%s %s %s %d\n", it->first.c_str(), _it->first.c_str(), _it->second.c_str(), net.minRoutingLayer+1);
	}
    }
    std::unordered_map<std::string, std::vector<std::pair<Point,Point>>>& routingSegments = routingGraph->get_routingSegments();
    for(auto it = routingSegments.begin(); it != routingSegments.end(); ++it) {
        for(auto _it = it->second.begin(); _it != it->second.end(); ++_it) {
	    printf("%s %d %d %d %d %d %d\n", it->first.c_str(), _it->first.x+1, _it->first.y+1, _it->first.z+1, _it->second.x+1, _it->second.y+1, _it->second.z+1);
	}
    }
}

void calculate_demand(RoutingGraph* routingGraph)
{
    int row = routingGraph->get_row();
    int column = routingGraph->get_column();
    int layer = routingGraph->get_layer();
    std::vector<std::vector<std::vector<Gcell>>>& grids = routingGraph->get_grids();
    std::unordered_map<std::string, Cell>& cellInstances = routingGraph->get_cellInstances();
    std::vector<std::vector<std::unordered_set<std::string>>>& placement = routingGraph->get_placement();
    std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<int, int>>>& sameGGrid = routingGraph->get_sameGGrid();
    std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<int, int>>>& adjHGGrid = routingGraph->get_adjHGGrid();
    std::unordered_map<std::string, std::vector<std::pair<Point,Point>>>& routingSegments = routingGraph->get_routingSegments();
    unordered_map<string, MasterCell> masterCells = routingGraph->get_masterCells();
    std::unordered_map<std::string, Net>& nets = routingGraph->get_nets();
    for(auto it = nets.begin(); it != nets.end(); ++it) {
        Net& net = it->second;
        for(auto _it = net.pins.begin(); _it != net.pins.end(); ++_it) {
            Cell& cell = cellInstances[_it->first];
	    int x = cell.x;
	    int y = cell.y;
	    int pinLayer = masterCells[cell.masterCellName].pins[_it->second].layer;
	    grids[x][y][pinLayer].insert_net(it->first);
        }
    }   
    for(int i = 0; i < row; i++) {
        for(int j = 0;j < column; j++) {
	    unordered_set<string> clCur = placement[i][j];
	    unordered_set<string> clPre, clNxt;
	    if(j > 0) clPre = placement[i][j-1];
	    if(j < column-1) clNxt = placement[i][j+1];
	    unordered_map<string,int> mcCntCur, mcCntPre, mcCntNxt;
	    for(auto it = clCur.begin(); it != clCur.end(); ++it) {
	        mcCntCur[cellInstances[*it].masterCellName]++;
		MasterCell& masterCell = masterCells[cellInstances[*it].masterCellName];
		for(auto _it = masterCell.blockages.begin(); _it != masterCell.blockages.end(); ++_it) {
		    grids[i][j][_it->second.layer].adjust_demand(_it->second.demand);
		}
	    }
	    for(auto it = clPre.begin(); it != clPre.end(); ++it) {
	        mcCntPre[cellInstances[*it].masterCellName]++;
            }
	    for(auto it = clNxt.begin(); it != clNxt.end(); ++it) {
	        mcCntNxt[cellInstances[*it].masterCellName]++;
	    }
	    for(auto it = sameGGrid.begin(); it != sameGGrid.end(); ++it) {
	        for(auto _it = it->second.begin(); _it != it->second.end(); ++_it) {
		    for(auto __it = _it->second.begin(); __it != _it->second.end(); ++__it) {
		        string MC1 = it->first, MC2 = _it->first;
			int layer = __it->first, extraDemand = __it->second;
			int pairCnt = min(mcCntCur[MC1], mcCntCur[MC2]);
			grids[i][j][layer].adjust_demand(extraDemand*pairCnt);
		    }
		}
	    }
	    for(auto it = adjHGGrid.begin(); it != adjHGGrid.end(); ++it) {
                for(auto _it = it->second.begin(); _it != it->second.end(); ++_it) {
                    for(auto __it = _it->second.begin(); __it != _it->second.end(); ++__it) {
			string MC1 = it->first, MC2 = _it->first;
                        int layer = __it->first, extraDemand = __it->second;
                        int pairCntPre = 0, pairCntNxt = 0;
			if(MC1 == MC2) {
			    pairCntPre = min(mcCntCur[MC1], mcCntPre[MC1]);
			    pairCntNxt = min(mcCntCur[MC1], mcCntNxt[MC1]);
			}
			else {
			    pairCntPre = min(mcCntCur[MC1], mcCntPre[MC2]) + min(mcCntCur[MC2], mcCntPre[MC1]);
			    pairCntNxt = min(mcCntCur[MC1], mcCntNxt[MC2]) + min(mcCntCur[MC2], mcCntNxt[MC1]);
			}
                        grids[i][j][layer].adjust_demand(extraDemand*(pairCntPre + pairCntNxt));
                    }
                }
            }
	}
    }
    for(auto it = routingSegments.begin(); it != routingSegments.end(); ++it) {
        for(auto _it = it->second.begin(); _it != it->second.end(); ++_it) {
	    int s_x = _it->first.x, s_y = _it->first.y, s_z = _it->first.z, e_x = _it->second.x, e_y = _it->second.y, e_z = _it->second.z;
	    if(s_x != e_x) {
	        for(int i = min(s_x,e_x); i <= max(s_x,e_x); i++) {
		    grids[i][s_y][s_z].insert_net(it->first);
		}
	    }
	    else if(s_y != e_y) {
	        for(int i = min(s_y,e_y); i <= max(s_y,e_y); i++) {
                    grids[s_x][i][s_z].insert_net(it->first);
                }
	    }
	    else if(s_z != e_z) {
	        for(int i = min(s_z,e_z); i <= max(s_z,e_z); i++) {
                    grids[s_x][s_y][i].insert_net(it->first);
                }
	    }
	}
    }

    for(int i = 0;i < row;i++) {
        for(int j = 0; j < column; j++) {
	    for(int k = 0; k < layer;k++) {
		grids[i][j][k].add_net_demand();
		
                printf("%d %d %d %d\n", i+1,j+1,k+1,grids[i][j][k].get_demand());
	    }
	}
    }
}

int main(int argc, char* argv[])
{
    RoutingGraph* routingGraph;
    string caseFile = argv[1];
    routingGraph = parser(routingGraph, caseFile);
    //debug_print(routingGraph);
    calculate_demand(routingGraph);
    return 0;
}


