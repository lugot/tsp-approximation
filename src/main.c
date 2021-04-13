#include "tsp.h"
#include "parsers.h"
#include "utils.h"
#include "solvers.h"
#include <stdio.h>
#include <cplex.h>
#include <stdlib.h>

int main(int argc, char** argv) {

    cplex_params params = (cplex_params) calloc(1, sizeof(struct cplex_params_t));
    run_options options = (run_options) calloc(1, sizeof(struct run_options_t));
    parse_command_line(argc, argv, params, options);

    printf("main: input file %s\n", options->model_name);

    instance inst = parse_input_file(options->model_name, "tsp", TSPLIB);
    add_params(inst, params);

    print_instance(inst, 1);

    solution sol = TSPopt(inst, SYMMETRIC_BENDERS_CALLBACK);

    print_solution(sol, 1);
    plot_solution_graphviz(sol);


    return EXIT_SUCCESS;
}
