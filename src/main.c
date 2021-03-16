#include "tsp.h"
#include "parsers.h"
#include "solvers.h"
#include <stdio.h>
#include <cplex.h>
#include <stdlib.h>

int main(int argc, char** argv) {

    instance inst = create_tsp_instance();
    print_instance(inst);

    parse_command_line(argc, argv, inst);
    parse_input(inst);
    print_instance(inst);

    solution optimal = TSPopt(inst);

    print_solution(optimal);
    plot_solution_graphviz(inst, optimal);

    return EXIT_SUCCESS;
}
