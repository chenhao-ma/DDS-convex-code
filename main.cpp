#include <iostream>
#include "Graph.h"
#include "Args.h"


int main(int argc, char *argv[]) {
    setbuf(stdout, NULL);
    clock_t begin = clock();

    Args *args = new Args();
    args->parse_args(argc, argv);
    FILE* d_file = fopen(args->address, "r");
    Graph g = Graph(d_file, args->NT, args->epsilon, args->divide_by_number, args->restricted, false);

    clock_t io_end = clock();

    printf("ok\n");
    if (args->exact) {
        printf("density %.4f\n", g.exact());
    } else if (args->stoc2020){
        printf("density %.4f\n", g.stoc2020());
    } else {
        printf("density %.4f\n", g.approx());
    }
    g.output_ds(args->ds_address);

    clock_t end = clock();
    double io_secs = double(io_end - begin) / CLOCKS_PER_SEC;
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    printf("io time: %.4f, total time: %.4f\n", io_secs, elapsed_secs);

    return 0;
}
