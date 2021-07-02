//
// Created by MA Chenhao on 11/11/2020.
//

#ifndef DDSAPP_EDGE_H
#define DDSAPP_EDGE_H

#include <utility>
using namespace std;

class Edge {
public:
    pair<int, int> nodes;
    pair<double, double> alpha;
    bool flag;
    Edge(int u, int v);
};


#endif //DDSAPP_EDGE_H
