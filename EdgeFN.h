//
// Created by MA Chenhao on 27/11/2020.
//

#ifndef DDSAPP_EDGEFN_H
#define DDSAPP_EDGEFN_H


class EdgeFN {
public:
    int from, to, index;
    double cap, flow;
    EdgeFN(int from, int to, double cap, double flow, double index);
};


#endif //DDSAPP_EDGEFN_H
