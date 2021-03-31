#include "tsp.h"
#include "parsers.h"
#include "utils.h"
#include "solvers.h"
#include <stdio.h>
#include <cplex.h>
#include <stdlib.h>

int main(int argc, char** argv) {

    instance inst = create_tsp_instance();
    /*print_instance(inst);*/

    parse_command_line(argc, argv, inst);
    parse_input_file(inst, "tsp");
    print_instance(inst, 1);

    TSPopt(inst, ASYMMETRIC_MTZ);

    print_instance(inst, 0);
    plot_solutions_graphviz(inst->sols, inst->num_solutions);
    /*plot_solution_graphviz(inst->sols[2]);*/


    return EXIT_SUCCESS;
}
