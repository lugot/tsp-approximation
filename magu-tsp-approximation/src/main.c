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
    parse_input_file(inst);
    print_instance(inst);

    solution optimal_cplex = TSPopt(inst);

    print_solution(optimal_cplex);
    solution tour = parse_optimal_tour(inst);
    /*plot_solutions_graphviz(inst, (solution []) {optimal_cplex, tour}, 2, 1);*/

    plot_solution_graphviz(inst, optimal_cplex);


    return EXIT_SUCCESS;
}
