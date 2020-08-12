#include<bits/stdc++.h>
#ifndef _RSMT_
#define _RSMT_

class RoutingGraph;

class RSMT{
public:
    RSMT() = delete;
    RSMT(RoutingGraph& g);
    RoutingGraph& routingGraph;
    void RSMT_reroute();
};

#endif
