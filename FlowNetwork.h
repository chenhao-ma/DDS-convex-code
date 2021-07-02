//
// Created by MA Chenhao on 27/11/2020.
//

#ifndef DDSAPP_FLOWNETWORK_H
#define DDSAPP_FLOWNETWORK_H

#include <vector>
#include <queue>
#include "EdgeFN.h"
#include <cmath>
#include <map>

using namespace std;

class FlowNetwork {
public:
    int n;
    vector <vector <EdgeFN>> adj;
    vector <double> excess;
    vector <int> dist, count;
    vector <bool> active;
    vector <vector <int>> B;
    vector<int> nums;
    map<int, int> mapping[2];
    int b;
    double m;

    FlowNetwork(vector<pair<int, int>> edges, double g, double sqrt_c);

    void add_edge(int from, int to, double cap);

    void enqueue (int v);

    void push (EdgeFN &e);

    void gap (int k);

    void relabel (int v);

    void discharge (int v);

    double get_maxflow(int s, int t, bool need_initial = true);
};


#endif //DDSAPP_FLOWNETWORK_H
