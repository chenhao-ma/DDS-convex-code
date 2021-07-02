#!/usr/bin/python3

def calc_f1(S, T, S_app, T_app):
    if len(S_app) + len(T_app):
        precision = (len(S & S_app) + len(T & T_app)) * 1.0 / (len(S_app) + len(T_app))
    else:
        precision = 0
    recall = (len(S & S_app) + len(T & T_app)) * 1.0 / (len(S) + len(T))
    if precision + recall > 0:
        f1 = 2 * precision * recall / (precision + recall)
    else:
        f1 = 0
    return f1, precision, recall


def read_S_T(address):
    S = set([])
    T = set([])
    try:
        file = open(address, "r")
        s, t = list(map(int, file.readline().strip().split(' ')))
        for i in range(s):
            S.add(int(file.readline().strip()))

        for i in range(t):
            T.add(int(file.readline().strip()))
    except OSError as e:
        print(e.errno)

    return S, T


# for dataset in ["MO", "TC", "OF", "AD", "AM"]:
for dataset in ["MO", "AD"]:
    exact_add = "output/" + dataset + "-e-0.0.txt"
    S, T = read_S_T(exact_add)
    for app in ["-a-0.1.txt", "-a-1.0.txt", "-a-2.0.txt"]:
        app_add = "output/" + dataset + app
        S_app, T_app = read_S_T(app_add)
        f1, precision, recall = calc_f1(S, T, S_app, T_app)
        print(dataset + app, f1, precision, recall)
