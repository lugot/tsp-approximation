#include "tsp.h"
#include "parser_utils.h"
#include <stdio.h>
#include <cplex.h>

int main(int argc, char **argv) {

    tsp_instance inst = create_tsp_instance();
    print_instance(inst);

    parse_command_line(argc, argv, inst);
    print_instance(inst);

    parse_input(inst);
    print_instance(inst);

    return 0;
}
