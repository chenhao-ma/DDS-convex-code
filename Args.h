//
// Created by MA Chenhao on 29/9/2020.
//

#ifndef DDSAPP_ARGS_H
#define DDSAPP_ARGS_H

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <getopt.h>

#define USAGE_TXT							   \
    "usage: \n"                                \
    "\t[-g input directed graph file]\n"       \
    "\t[-t repeated iterations]\n"             \
    "\t[-e epsilon]\n"                         \
    "\t[-a approximate or exact DDS or stoc2020]\n"        \
    "\t[-d n by number, v by value]\n"

class Args {
public:
    char *address{};
    char ds_address[50];
    int NT = 100;
    double epsilon = 1;
    bool exact = false;
    bool divide_by_number = true;
    bool stoc2020 = false;
    bool restricted = false;
    ~Args();
    void usage(char *msg, int exit_status);
    void parse_args(int argc, char *argv[]);
};


#endif //DDSAPP_ARGS_H
