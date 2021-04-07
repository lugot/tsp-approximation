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

<<<<<<< HEAD
    double runtime = 0.0;

    for (int i=0; i<5; i++) {
        runtime += TSPopt(inst, SYMMETRIC_BENDERS)->solve_time;
        /*print_instance(inst, 1);*/
        printf("heo\n");
    }

    printf("time: %lf\n", runtime/1000.0);
    /*print_solution(inst->sols[0], 1);*/
    /*plot_solutions_graphviz(inst->sols, inst->num_solutions);*/
    /*plot_solution_graphviz(inst->sols[0]);*/
=======
    TSPopt(inst, ASYMMETRIC_MTZ);

    print_instance(inst, 0);
    plot_solutions_graphviz(inst->sols, inst->num_solutions);
    /*plot_solution_graphviz(inst->sols[2]);*/

>>>>>>> parent of f2d3960 (added benders)

    return EXIT_SUCCESS;
}
