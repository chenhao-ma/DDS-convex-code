//
// Created by MA Chenhao on 11/11/2020.
//

#include "Edge.h"


Edge::Edge(int u, int v) {
    nodes.first = u;
    nodes.second = v;
    alpha.first = alpha.second = 0;
    flag = true;
}
