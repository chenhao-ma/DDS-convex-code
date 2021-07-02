//
// Created by MA Chenhao on 28/9/2020.
//

#ifndef DDSAPP_GRAPH_H
#define DDSAPP_GRAPH_H

#include <vector>
#include <cstdio>
#include <tuple>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <set>
#include "Edge.h"
#include "FlowNetwork.h"

using namespace std;

class Graph {
public:
    Graph(FILE *pFile, int NT, double epsilon, bool by_number, bool res, bool ablation);

    double approx();
    double exact();
    double stoc2020();
    void output_ds(char *ds_address);

private:
    int n;
    int m;
    int NT;
    double epsilon;
    bool divide_by_number = true;
    bool restricted = false;
    double res_width = 2;
    bool ablation_test = false;
    int CT;
    vector<vector<vector<int>>> adj;
    vector<vector<double>> r;
    vector<vector<int>> deg;
    vector<vector<bool>> selected;
    vector<int> slt_edges;
    vector<Edge> edges;
    vector<pair<int, int>> slt_nodes;
    vector<int> S;
    vector<int> T;
    double density;
    vector<int> max_deg;
    vector<int> max_pos;

    void frank_wolfe(int CT, double c);
    tuple<double, double, bool> app_cDDS(double c);
    tuple<double, double, bool> exact_cDDS(double c);
    bool stoc_vwDS(double sqrt_c);
    void divide_conquer(double c_l, double c_r);
    void can_nodes(double c_l, double c_r);
    void can_nodes_core(double c_l, double c_r);
    void init_DDS();
};


#endif //DDSAPP_GRAPH_H
