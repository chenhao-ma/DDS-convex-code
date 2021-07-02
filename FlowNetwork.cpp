//
// Created by MA Chenhao on 27/11/2020.
//

#include "FlowNetwork.h"

FlowNetwork::FlowNetwork(vector<pair<int, int>> edges, double g, double sqrt_c) {
    n = 0;
    for (auto & edge : edges) {
        if (mapping[0].find(edge.first) == mapping[0].end()) {
            ++n;
            mapping[0].insert(make_pair(edge.first, n));
        }
    }
    nums.push_back(n);
    for (auto & edge : edges) {
        if (mapping[1].find(edge.second) == mapping[1].end()) {
            ++n;
            mapping[1].insert(make_pair(edge.second, n));
        }
    }
    nums.push_back(n);
    n += 2;
    adj.resize(n);

    for (auto &edge : edges) {
        int u = mapping[1].find(edge.second)->second;
        int v = mapping[0].find(edge.first)->second;
        add_edge(u, v, 2);
    }

    for (int u = nums[0] + 1; u <= nums[1]; u++) {
        add_edge(0, u, adj[u].size());
    }

    for (int u = 1; u < n - 1; u++) {
        if (u <= nums[0]) {
            add_edge(u, n - 1, g / 2 / sqrt_c);
        } else {
            add_edge(u, n - 1, g * sqrt_c / 2);
        }
    }
}

void FlowNetwork::add_edge(int from, int to, double cap) {
    adj[from].push_back(EdgeFN(from, to, cap, 0, adj[to].size()));
    if (from == to) {
        adj[from].back().index++;
    }
    adj[to].push_back(EdgeFN(to, from, 0, 0, adj[from].size() - 1));
}

void FlowNetwork::enqueue (int v) {
    if (!active[v] && excess[v] > 0 && dist[v] < n) {
        active[v] = true;
        B[dist[v]].push_back(v);
        b = std::max(b, dist[v]);
    }
}

void FlowNetwork::push (EdgeFN &e) {
    double amt = std::min(excess[e.from], e.cap - e.flow);
    if (dist[e.from] == dist[e.to] + 1 && amt > 0) {
        e.flow += amt;
        adj[e.to][e.index].flow -= amt;
        excess[e.to] += amt;
        excess[e.from] -= amt;
        enqueue(e.to);
    }
}

void FlowNetwork::gap (int k) {
    for (int v = 0; v < n; v++) {
        if (dist[v] >= k) {
            count[dist[v]]--;
            dist[v] = std::max(dist[v], n);
            count[dist[v]]++;
            enqueue(v);
        }
    }
}

void FlowNetwork::relabel(int v) {
    count[dist[v]]--;
    dist[v] = n;
    for (auto e : adj[v]) if (e.cap - e.flow > 0) {
            dist[v] = std::min(dist[v], dist[e.to] + 1);
        }
    count[dist[v]]++;
    enqueue(v);
}

void FlowNetwork::discharge(int v) {
    for (auto &e : adj[v]) {
        if (excess[v] > 0) {
            push(e);
        } else {
            break;
        }
    }

    if (excess[v] > 0) {
        if (count[dist[v]] == 1) {
            gap(dist[v]);
        } else {
            relabel(v);
        }
    }
}

double FlowNetwork::get_maxflow(int s, int t, bool need_initial) {
    if (need_initial) {
        dist = std::vector<int>(n, 0);
        excess = std::vector<double>(n, 0);
        count = std::vector<int>(n + 1, 0);
        active = std::vector<bool>(n, false);
        B = std::vector<std::vector<int>>(n);
        b = 0;

        for (auto &e: adj[s]) {
            excess[s] += e.cap;
        }

        count[0] = n;
        enqueue(s);
        active[t] = true;
    }

    while (b >= 0) {
        if (!B[b].empty()) {
            int v = B[b].back();
            B[b].pop_back();
            active[v] = false;
            discharge(v);
        } else {
            b--;
        }
    }

    return excess[t];
}
