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
    printf("here\n");
    /*parse_input_file(inst, "tsp");*/
    parse_input_file(inst, "tsp");
    print_instance(inst);
    /*print_instance(inst);*/

    solution optimal_cplex = TSPopt(inst, ASYMMETRIC_MTZ);

    print_solution(optimal_cplex);
    plot_solution_graphviz(inst->sols[0]);


    return EXIT_SUCCESS;
}
