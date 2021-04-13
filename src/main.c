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

    solution sol = TSPopt(inst, ASYMMETRIC_GG);

    print_solution(sol, 1);
    plot_solution_graphviz(sol);

    /*instance* a = create_random_instances(10, 20, 3);*/
    /*for (int i=0; i<10; i++) print_instance(a[i], 1);*/

    /*instance inst = create_tsp_instance();*/
    /*[>print_instance(inst);<]*/

    /*parse_command_line(argc, argv, inst);*/
    /*parse_input_file(inst, "tsp");*/
    /*[>print_instance(inst, 1);<]*/

    /*TSPopt(inst, SYMMETRIC_BENDERS);*/

    /*print_instance(inst, 1);*/
    /*[>print_solution(inst->sols[0], 1);<]*/
    /*[>plot_solutions_graphviz(inst->sols, inst->num_solutions);<]*/
    /*plot_solution_graphviz(inst->sols[0]);*/

    return EXIT_SUCCESS;
}
