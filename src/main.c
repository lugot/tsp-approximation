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

    instance dummy = duplicate_instance_parameters(inst);
    parse_input_file(dummy, "opt.tour");

    TSPopt(inst, ASYMMETRIC_MTZ);
    zstar(inst, dummy->sols[0]);
    add_solution(inst, dummy->sols[0]);

    print_instance(inst, 0);
    plot_solution_graphviz(inst->sols[0]);


    return EXIT_SUCCESS;
}
