//
// Created by MA Chenhao on 29/9/2020.
//

#include "Args.h"

Args::~Args() {
    delete address;
}

void Args::usage(char *msg, int exit_status) {
    fprintf(exit_status == 0 ? stdout : stderr, "%s", USAGE_TXT);

    if (msg) {
        fprintf(exit_status == 0 ? stdout : stderr, "\n%s\n", msg);
    }
    exit(exit_status);
}

void Args::parse_args(int argc, char **argv) {
    int c;
    opterr = 0;
    if (argc < 2) {
        usage(nullptr, 0);
    }

    while ((c = getopt(argc, argv, "g:t:e:a:d:r")) != -1) {
        switch (c) {
            case 'g': {
                printf("%s\n", optarg);
                address = strdup(optarg);
                break;
            }

            case 't': {
                sscanf(optarg, "%d", &NT);
                break;
            }

            case 'e': {
                sscanf(optarg, "%lf", &epsilon);
                break;
            }

            case 'a':{
                char* tmp;
                tmp = strdup(optarg);
                exact = tmp[0] == 'e';
                stoc2020 = tmp[0] == 's';
                delete tmp;
                break;
            }

            case 'd':{
                char* tmp;
                tmp = strdup(optarg);
                divide_by_number = tmp[0] == 'n';
                delete tmp;
                break;
            }

            case 'r':{
                restricted = true;
                break;
            }

            default:
                break;
        }
    }
    char dataset[5];
    strncat(dataset, address, 2);
    if (exact) {
        sprintf(ds_address, "output/%s-e-%.1f.txt", dataset, epsilon);
    } else if (stoc2020) {
        sprintf(ds_address, "output/%s-s-%.1f.txt", dataset, epsilon);
    } else {
        sprintf(ds_address, "output/%s-a-%.1f.txt", dataset, epsilon);
    }

    printf("%s\n", ds_address);
}