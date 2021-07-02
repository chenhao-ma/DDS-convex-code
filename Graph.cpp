//
// Created by MA Chenhao on 28/9/2020.
//



#include "Graph.h"

Graph::Graph(FILE *file, int numT, double epsilon, bool by_number, bool res, bool ablation) {
    //    clock_t begin = clock();
    fscanf(file, "%d%d", &n, &m);
    divide_by_number = by_number;
    restricted = res;
    ablation_test = ablation;
    res_width = max(epsilon, res_width);

    //initialize adjacent list
    adj.resize(2);
    r.resize(2);
    deg.resize(2);
    selected.resize(2);
    for (int i = 0; i < 2; i++) {
        adj[i].resize(n);
        r[i].resize(n);
        deg[i].resize(n);
        selected[i].resize(n);
    }

    //read edges from the file
    for (int i = 0; i < m; i++) {
        int u, v;
        fscanf(file, "%d%d", &u, &v);
        edges.emplace_back(u, v);
        adj[0][u].push_back(edges.size() - 1);
        adj[1][v].push_back(edges.size() - 1);
    }

    max_deg.resize(2);
    max_pos.resize(2);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < n; j++) {
            if (adj[i][j].size() > max_deg[i]) {
                max_deg[i] = (int) adj[i][j].size();
                max_pos[i] = j;
            }
        }
    }

    //initialize maximal density
    density = 0;
    NT = numT;
    this->epsilon = epsilon;
}

void Graph::output_ds(char *ds_address) {
//    printf("S:\n");

    FILE* dsFile = fopen(ds_address, "w");
    fprintf(dsFile, "%lu %lu\n", S.size(), T.size());
//    fprintf(dsFile, "S:\n");
    for (auto u : S){
        fprintf(dsFile, "%d\n", u);
    }
//    fprintf(dsFile, "T:\n");
    for (auto u : T){
        fprintf(dsFile, "%d\n", u);
    }

    fclose(dsFile);
}

void Graph::frank_wolfe(int CT, double c) {
    if (CT == 0) {
        for (auto & edge : slt_edges) {
            edges[edge].alpha.first = edges[edge].alpha.second = 0.5;
        }
        for (auto & node : slt_nodes) {
            r[node.first][node.second] = 0;
            for (auto & i : adj[node.first][node.second]) {
                if (edges[i].flag) {
                    if (node.first == 0) {
                        r[node.first][node.second] += 2 * sqrt(c) * edges[i].alpha.first;
                    } else {
                        r[node.first][node.second] += 2 / sqrt(c) * edges[i].alpha.second;
                    }
                }
            }
        }
//        for (int i = 0; i < n; i++) {
//            r[0][i] = 0;
//            for (auto & j : adj[0][i]) {
//                r[0][i] += 2 * sqrt(c) * edges[j].alpha.first;
//            }
//        }
//        for (int i = 0; i < n; i++) {
//            r[1][i] = 0;
//            for (auto & j : adj[1][i]) {
//                r[1][i] += 2 / sqrt(c) * edges[j].alpha.second;
//            }
//        }
    }

    for (int t = CT + 1; t <= CT + NT; t++) {
        double gamma_t = 2.0 / (t + 2);
        for (auto & edge : slt_edges) {
            edges[edge].alpha.first *= 1 - gamma_t;
            edges[edge].alpha.second *= 1 - gamma_t;
            if (r[0][edges[edge].nodes.first] < r[1][edges[edge].nodes.second]) {
                edges[edge].alpha.first += gamma_t;
            } else {
                edges[edge].alpha.second += gamma_t;
            }
        }

        for (auto & node : slt_nodes) {
            r[node.first][node.second] = 0;
            for (auto & i : adj[node.first][node.second]) {
                if (edges[i].flag) {
                    if (node.first == 0) {
                        r[node.first][node.second] += 2 * sqrt(c) * edges[i].alpha.first;
                    } else {
                        r[node.first][node.second] += 2 / sqrt(c) * edges[i].alpha.second;
                    }
                }
            }
        }
    }
}

tuple<double, double, bool> Graph::app_cDDS(double c) {
//    slt_nodes.clear();
//    for (int i = 0; i < 2; i++) {
//        for (int j = 0; j < n; j++) {
//            slt_nodes.emplace_back(i, j);
//        }
//    }

    sort(slt_nodes.begin(), slt_nodes.end(), [this](pair<int, int> a, pair<int, int> b)->bool {
        return r[a.first][a.second] > r[b.first][b.second];
    });

//    if (r[slt_nodes[0].first][slt_nodes[0].second] < density * sqrt(1 + epsilon)) {
//        double t = r[slt_nodes[0].first][slt_nodes[0].second] / density / sqrt(1 + epsilon);
    double adjust = epsilon;
    if (adjust > 0.5) adjust = 0;
    if (r[slt_nodes[0].first][slt_nodes[0].second] < density * sqrt(1 + adjust)) {
        double t = r[slt_nodes[0].first][slt_nodes[0].second] / density / sqrt(1 + adjust);
        //x -> (2 c - c t^2 - 2 Sqrt[c^2 - c^2 t^2])/t^2, x -> (2 c - c t^2 + 2 Sqrt[c^2 - c^2 t^2])/t^2
        double c_o = (2 * c - c * t * t - 2 * sqrt(c * c - c * c * t * t)) / (t * t);
        double c_p = c * c / c_o;
        return make_tuple(c_o, c_p, true);
    }

    bool flag = false;
//    set<int> Sc, Tc;
    int edge_num = 0, best_pos = 0, pos = 0;
    double rho_c = 0, rho = 0, best_c_prime = 0;
    for (auto & node : slt_nodes) {
        selected[node.first][node.second] = false;
    }
    vector<int> cnt;
    cnt.resize(2, 0);
    for (auto & node : slt_nodes) {
        ++pos;
        if (node.first == 0) {
            ++cnt[0];
//            Sc.insert(node.second);
        } else {
            ++cnt[1];
//            Tc.insert(node.second);
        }
        selected[node.first][node.second] = true;
//        if (Sc.empty() || Tc.empty()) continue;
        if (cnt[0] == 0 || cnt[1] == 0) continue;
        for (auto & i : adj[node.first][node.second]) {
//            if ((node.first == 0 ? Tc.count(edges[i].nodes.second) : Sc.count(edges[i].nodes.first)) > 0) ++edge_num;
            if (node.first == 0 ? selected[1][edges[i].nodes.second] : selected[0][edges[i].nodes.first]) ++edge_num;
        }

//        if (edge_num / sqrt(Sc.size() * Tc.size()) > rho) {
//            rho = edge_num / sqrt(Sc.size() * Tc.size());
//            best_pos = pos;
//        }
        if (edge_num / sqrt( (double) cnt[0] * cnt[1]) > rho) {
            rho = edge_num / sqrt( (double) cnt[0] * cnt[1]);
            best_pos = pos;
//            printf("%d %d %d %.4f\n", edge_num, cnt[0], cnt[1], rho);
        }
//        double c_prime = (double) Sc.size() / Tc.size();
        double c_prime = (double) cnt[0] / cnt[1];
//        rho_c = 2 * sqrt(c * c_prime) / (c + c_prime) * edge_num / sqrt(Sc.size() * Tc.size());
        rho_c = 2 * sqrt(c * c_prime) / (c + c_prime) * edge_num / sqrt( (double) cnt[0] * cnt[1]);
        if (r[slt_nodes[0].first][slt_nodes[0].second] / rho_c < 1 + epsilon) {
            if (!flag || abs(c_prime - c) > abs(best_c_prime - c)) {
                best_c_prime = c_prime;
            }
            flag = true;
        }
//        double c_prime = (double) Sc.size() / Tc.size();
//        if (2 * sqrt(c * c_prime) / (c + c_prime) * edge_num / sqrt(Sc.size() * Tc.size()) > rho_c) {
//            best_c_prime = c_prime;
//            rho = edge_num / sqrt(Sc.size() * Tc.size());
//            rho_c = 2 * sqrt(c * c_prime) / (c + c_prime) * rho;
//            best_pos = pos;
//        }
    }

    if (rho > density) {
        density = rho;
        S.clear(); T.clear();
        for (int i = 0; i < best_pos; i++) {
            if (slt_nodes[i].first == 0)
                S.push_back(slt_nodes[i].second);
            else
                T.push_back(slt_nodes[i].second);
        }
    }

    double c_o = best_c_prime;
    double c_p = c * c / c_o;
    if (!flag) return make_tuple(c_o, c_p, flag);
    if (c_o > c_p) swap(c_o, c_p);
    double cl = c / (1 + epsilon), cr = c * (1 + epsilon);
    if (c_o <= cl && c_p >= cr) return make_tuple(c_o, c_p, flag);
    if (r[slt_nodes[0].first][slt_nodes[0].second] / rho_c < sqrt(1 + epsilon))
        return make_tuple(min(c_o, cl), max(c_p, cr), flag);
    else
        return make_tuple(c_o, c_p, false);

//    return make_tuple(c_o, c_p, r[slt_nodes[0].first][slt_nodes[0].second] / rho_c < 1 + epsilon);
}

tuple<double, double, bool> Graph::exact_cDDS(double c) {
    //sorting the nodes according to the descending order of r
    sort(slt_nodes.begin(), slt_nodes.end(), [this](pair<int, int> a, pair<int, int> b)->bool {
        return r[a.first][a.second] > r[b.first][b.second];
    });

//    printf("first r value %.4f, density %.4f\n", r[slt_nodes[0].first][slt_nodes[0].second], density);
    if (r[slt_nodes[0].first][slt_nodes[0].second] < density) {
        double t = r[slt_nodes[0].first][slt_nodes[0].second] / density;
        //x -> (2 c - c t^2 - 2 Sqrt[c^2 - c^2 t^2])/t^2, x -> (2 c - c t^2 + 2 Sqrt[c^2 - c^2 t^2])/t^2
        double c_o = (2 * c - c * t * t - 2 * sqrt(c * c - c * c * t * t)) / (t * t);
        double c_p = c * c / c_o;
        return make_tuple(c_o, c_p, true);
    }

    //tentative cDDS extracting
//    set<int> Sc, Tc;
    vector<int> cnt; cnt.resize(2, 0);
    int edge_num = 0, best_pos = 0, pos = 0;
    double rho_c = 0, rho = 0, best_c_prime = 0;
    for (auto & node : slt_nodes) {
        selected[node.first][node.second] = false;
    }
    for (auto & node : slt_nodes) {
        ++pos;
        if (node.first == 0) {
//            Sc.insert(node.second);
            ++cnt[0];
        } else {
            ++cnt[1];
//            Tc.insert(node.second);
        }
        selected[node.first][node.second] = true;
//        if (Sc.empty() || Tc.empty()) continue;
        if (cnt[0] == 0 || cnt[1] == 0) continue;
        for (auto & i : adj[node.first][node.second]) {
//            if ((node.first == 0 ? Tc.count(edges[i].nodes.second) : Sc.count(edges[i].nodes.first)) > 0) ++edge_num;
            if (node.first == 0 ? selected[1][edges[i].nodes.second] : selected[0][edges[i].nodes.first]) ++edge_num;
        }
//        double c_prime = (double) Sc.size() / Tc.size();
        double c_prime = (double) cnt[0] / cnt[1];
//        if (2 * sqrt(c * c_prime) / (c + c_prime) * edge_num / sqrt(Sc.size() * Tc.size()) > rho_c) {
        if (2 * sqrt(c * c_prime) / (c + c_prime) * edge_num / sqrt((double) cnt[0] * cnt[1]) > rho_c) {
            best_c_prime = c_prime;
//            rho = edge_num / sqrt(Sc.size() * Tc.size());
            rho = edge_num / sqrt((double) cnt[0] * cnt[1]);
            rho_c = 2 * sqrt(c * c_prime) / (c + c_prime) * rho;
            best_pos = pos;
        }
    }

    printf("r-infty %.4f rho-c %.4f\n", r[slt_nodes[0].first][slt_nodes[0].second], rho_c);

    bool flag = true;
    if (ablation_test) {
        if (abs(r[slt_nodes[0].first][slt_nodes[0].second] - rho_c) <= 1e-3) {
            flag = true;
        } else {
            flag = false;
        }
    } else {
//Stable testing
//    Sc.clear(); Tc.clear();
//    for (int i = 0; i < best_pos; i++) {
//        if (slt_nodes[i].first == 0)
//            Sc.insert(slt_nodes[i].second);
//        else
//            Tc.insert(slt_nodes[i].second);
//    }
        for (int i = best_pos; i < slt_nodes.size(); i++) {
            selected[slt_nodes[i].first][slt_nodes[i].second] = false;
        }

//    printf("S %lu T %lu\n", Sc.size(), Tc.size());

        for (auto & e : slt_edges) {
//        if (Sc.count(edges[e].nodes.first) > 0 && Tc.count(edges[e].nodes.second) == 0) {
            if (selected[0][edges[e].nodes.first] && !selected[1][edges[e].nodes.second]) {
                r[0][edges[e].nodes.first] -= 2.0 * sqrt(c) * edges[e].alpha.first;
                r[1][edges[e].nodes.second] += 2.0 / sqrt(c) * edges[e].alpha.first;
                edges[e].alpha.first = 0; edges[e].alpha.second = 1;
//        } else if (Sc.count(edges[e].nodes.first) == 0 && Tc.count(edges[e].nodes.second) > 0) {
            } else if (!selected[0][edges[e].nodes.first] && selected[1][edges[e].nodes.second]) {
                r[0][edges[e].nodes.first] += 2.0 * sqrt(c) * edges[e].alpha.second;
                r[1][edges[e].nodes.second] -= 2.0 / sqrt(c) * edges[e].alpha.second;
                edges[e].alpha.first = 1; edges[e].alpha.second = 0;
            }
        }

        double min_r = r[slt_nodes[0].first][slt_nodes[0].second];
        for (int i = 1; i < best_pos; ++i) {
            min_r = min(min_r, r[slt_nodes[i].first][slt_nodes[i].second]);
        }
        for (int i = best_pos; i < slt_nodes.size(); ++i) {
            if (min_r < r[slt_nodes[i].first][slt_nodes[i].second]) {
                double c_o = best_c_prime;
                double c_p = c * c / c_o;
                return make_tuple(c_o, c_p, false);
            }
        }

        //cDDS checking via edge num
        slt_nodes.resize(best_pos);

        //selected edges updating
        vector<int> new_slt_edges;
        for (auto & e : slt_edges) {
//        if (Sc.count(edges[e].nodes.first) > 0 && Tc.count(edges[e].nodes.second) > 0){
            if (selected[0][edges[e].nodes.first] && selected[1][edges[e].nodes.second]) {
                edges[e].flag = true;
                new_slt_edges.push_back(e);
            } else {
                edges[e].flag = false;
            }
        }
        slt_edges = new_slt_edges;

//    printf("CT %d #edge %lu stable \n", CT, slt_edges.size());

        //cDDS checking via max-flow
        vector<pair<int, int>> fn_edges;
        for (auto &e : slt_edges) {
            fn_edges.push_back(edges[e].nodes);
        }
        FlowNetwork fn = FlowNetwork(fn_edges, rho_c, sqrt(c));
        double max_flow = fn.get_maxflow(0, fn.n-1);
//    printf("rho %.6f, rho_c %.6f, max flow %.6f, slt edges %lu, delta %.6f\n", rho, rho_c, max_flow, slt_edges.size(), abs(max_flow - slt_edges.size()));
        if (abs(max_flow - slt_edges.size()) <= 1e-3) {
            flag = true;
        } else {
            flag = false;
        }
    }

    //updating DDS if possible
    if (rho > density) {
        density = rho;
        S.clear(); T.clear();
        for (int i = 0; i < best_pos; i++) {
            if (slt_nodes[i].first == 0)
                S.push_back(slt_nodes[i].second);
            else
                T.push_back(slt_nodes[i].second);
        }
    }
//    if (flag) {
//        printf("r-infty %.4f rho-c %.4f\n", r[slt_nodes[0].first][slt_nodes[0].second], rho_c);
//    }
    double c_o = best_c_prime;
    double c_p = c * c / c_o;
    return make_tuple(c_o, c_p, flag);
}

double Graph::approx() {
    init_DDS();
    divide_conquer(1.0/n, n);
    return density;
}

void Graph::divide_conquer(double c_l, double c_r) {
//    printf("init cl %.4f cr %.4f\n", c_l, c_r);
    c_l = max(c_l, (density / 2 / max_deg[0]) * (density / 2 / max_deg[0]));
    c_r = min(c_r, (2 * max_deg[1] / density) * (2 * max_deg[1] / density));
    if (c_l > c_r) return;
    double c = (c_l + c_r) / 2;
    if (divide_by_number) {
        if (c_l < 1 && c_r > 1) {
            c = 1;
        } else if (c_r <= 1) {
            c = (c_l + c_r) / 2;
        } else if (c_l >= 1) {
            c = 2 / (1 / c_l + 1 / c_r);
        }
    }
    printf("cl %.4f c %.4f cr %.4f\n", c_l, c, c_r);
    clock_t fw_start = clock();
    CT = 0;
    double c_o = 0, c_p = 0;
//    can_nodes(c_l, c_r);
    double res_c_l = c_l, res_c_r = c_r;
    if (restricted) {
        res_c_l = max(c_l, c / res_width);
        res_c_r = min(c_r, c * res_width);
    }
    //todo
    can_nodes_core(res_c_l, res_c_r);
    printf("#edge %lu\n", slt_edges.size());
    if (!slt_nodes.empty()) {
        while (true) {
            frank_wolfe(CT, c);
            auto ret = (epsilon > 0) ? app_cDDS(c) : exact_cDDS(c);
//        printf("%.4f %.4f %.4f\n", c, ret.first, ret.second);
            CT += NT;
            if (get<2>(ret)) {
                c_o = get<0>(ret);
                c_p = get<1>(ret);
                break;
            }
        }
    } else if (restricted) {
        c_o = res_c_l;
        c_p = res_c_r;
    } else {
        return;
    }
    printf("CT %d\n", CT);
    clock_t fw_end = clock();
    double fw_secs = double(fw_end - fw_start) / CLOCKS_PER_SEC;
    printf("fw time: %.4f\n", fw_secs);

    if (c_o > c_p) swap(c_o, c_p);
    if (restricted) {
        c_o = max(res_c_l, c_o);
        c_p = min(res_c_r, c_p);
    }
    if (c_l < c_o) divide_conquer(c_l, c_o);
    if (c_p < c_r) divide_conquer(c_p, c_r);
}

void Graph::can_nodes(double c_l, double c_r) {
    slt_nodes.clear();
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < n; j++) {
            if (i == 0){
                if (2 * sqrt(c_r) * adj[i][j].size() > density) {
                    slt_nodes.emplace_back(i, j);
                    selected[i][j] = true;
                } else {
                    selected[i][j] = false;
                }
            } else {
                if (2 / sqrt(c_l) * adj[i][j].size() > density) {
                    slt_nodes.emplace_back(i, j);
                    selected[i][j] = true;
                } else {
                    selected[i][j] = false;
                }
            }
        }
    }
    slt_edges.clear();
    for (int i = 0; i < m; i++) {
        if (selected[0][edges[i].nodes.first] && selected[1][edges[i].nodes.second]) {
            slt_edges.push_back(i);
            edges[i].flag = true;
        } else {
            edges[i].flag = false;
        }
    }
}

void Graph::can_nodes_core(double c_l, double c_r) {
    printf("updated cl %.4f cr %.4f\n", c_l, c_r);
    slt_nodes.clear();
    queue<pair<int, int> > q;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < n; j++) {
            deg[i][j] = adj[i][j].size();
            if (i == 0 && deg[i][j] * 2 * sqrt(c_r) <= density ||
                i == 1 && deg[i][j] * 2 / sqrt(c_l) <= density ) {
                q.push(make_pair(i, j));
                selected[i][j] = false;
            } else {
                selected[i][j] = true;
            }
//            if (i == 0 && deg[i][j] * 2 * sqrt(c_r) <= density * (1 + epsilon) ||
//                i == 1 && deg[i][j] * 2 / sqrt(c_l) <= density * (1 + epsilon)) {
//                q.push(make_pair(i, j));
//                selected[i][j] = false;
//            } else {
//                selected[i][j] = true;
//            }
        }
    }
    for (auto & edge : edges) {
        edge.flag = true;
    }
    while (!q.empty()) {
        auto node = q.front();
        q.pop();
        for (auto & e : adj[node.first][node.second]) {
            edges[e].flag = false;
            if (node.first == 0) {
                if (selected[1][edges[e].nodes.second]) {
                    --deg[1][edges[e].nodes.second];
                    if (deg[1][edges[e].nodes.second] * 2 / sqrt(c_l) <= density * (1 + epsilon)) {
                        selected[1][edges[e].nodes.second] = false;
                        q.push(make_pair(1, edges[e].nodes.second));
                    }
                }
            } else {
                if (selected[0][edges[e].nodes.first]) {
                    --deg[0][edges[e].nodes.first];
                    if (deg[0][edges[e].nodes.first] * 2 * sqrt(c_r) <= density * (1 + epsilon)) {
                        selected[0][edges[e].nodes.first] = false;
                        q.push(make_pair(0, edges[e].nodes.first));
                    }
                }
            }
        }
    }
    slt_edges.clear();
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < n; j++) {
            if (selected[i][j]) {
                slt_nodes.emplace_back(i, j);
                if (i == 0) {
                    for (auto & e : adj[i][j]) {
                        if (edges[e].flag) {
                            slt_edges.push_back(e);
                        }
                    }
                }
            }
        }
    }
}

double Graph::exact() {
    epsilon = 0;
    init_DDS();
    divide_conquer(1.0/n, n);
    return density;
}

void Graph::init_DDS() {
    int cur = (max_deg[0] > max_deg[1]) ? 0 : 1;
    density = sqrt(max_deg[cur]);
    S.clear(); T.clear();
    if (cur == 0) {
        S.push_back(max_pos[cur]);
        for (auto & e : adj[cur][max_pos[cur]]) {
            T.push_back(edges[e].nodes.second);
        }
    } else {
        T.push_back(max_pos[cur]);
        for (auto & e : adj[cur][max_pos[cur]]) {
            S.push_back(edges[e].nodes.first);
        }
    }

}

double Graph::stoc2020() {
    double eps = 1 - 1 / (1 + epsilon);
    printf("stoc2020 eps: %.4f\n", eps);
    double sqrt_c = 1 / sqrt(n);
    while (sqrt_c <= sqrt(n)) {
//        printf("sqrt_c %.4f\n", sqrt_c);
        slt_edges.clear();
        slt_nodes.clear();
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < n; j++) {
                selected[i][j] = true;
                slt_nodes.emplace_back(i, j);
            }
        }
        for (int i = 0; i < m; i++) {
            edges[i].flag = true;
            slt_edges.push_back(i);
        }
        int CT = 0;
        while (true) {
            frank_wolfe(CT, sqrt_c*sqrt_c);
            auto ret = stoc_vwDS(sqrt_c);
            if (ret) {
                break;
            }
            CT += NT;
        }
        sqrt_c *= 1 / (1 - eps / 2);
    }
    return density;
}

bool Graph::stoc_vwDS(double sqrt_c) {
    sort(slt_nodes.begin(), slt_nodes.end(), [this](pair<int, int> a, pair<int, int> b)->bool {
        return r[a.first][a.second] > r[b.first][b.second];
    });

    for (auto & node : slt_nodes) {
        selected[node.first][node.second] = false;
    }

    vector<int> cnt;
    cnt.resize(2, 0);
    double vw_rho = 0, rho = 0;
    int edge_num = 0, best_pos = 0, pos = 0;
    for (auto & node : slt_nodes) {
        ++pos;
        ++cnt[node.first];
        selected[node.first][node.second] = true;
        if (cnt[0] == 0 || cnt[1] == 0) {
            continue;
        }
        for (auto & i : adj[node.first][node.second]) {
//            if ((node.first == 0 ? Tc.count(edges[i].nodes.second) : Sc.count(edges[i].nodes.first)) > 0) ++edge_num;
            if (node.first == 0 ? selected[1][edges[i].nodes.second] : selected[0][edges[i].nodes.first])
                ++edge_num;
        }
        if (2 * edge_num / (cnt[0] / sqrt_c + cnt[1] * sqrt_c) > vw_rho) {
            vw_rho = 2 * edge_num / (cnt[0] / sqrt_c + cnt[1] * sqrt_c);
            rho = edge_num / sqrt( (double) cnt[0] * cnt[1]);
            best_pos = pos;
        }
    }
//    printf("r[0] %.4f, vw_rho %.4f, rho %.4f\n", r[slt_nodes[0].first][slt_nodes[0].second], vw_rho, rho);

    if (rho > density) {
        density = rho;
        S.clear(); T.clear();
        for (int i = 0; i < best_pos; i++) {
            if (slt_nodes[i].first == 0)
                S.push_back(slt_nodes[i].second);
            else
                T.push_back(slt_nodes[i].second);
        }
    }

    return vw_rho / r[slt_nodes[0].first][slt_nodes[0].second] >= 1 - epsilon / 2;
}


