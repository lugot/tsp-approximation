#include <assert.h>
#include <cplex.h>
#include <stdio.h>
#include <stdlib.h>

#include "../include/parsers.h"
#include "../include/solvers.h"
#include "../include/tsp.h"
#include "../include/utils.h"

int main(int argc, char** argv) {
    cplex_params params = create_params();
    run_options options = create_options();
    parse_command_line(argc, argv, params, options);

    if (options->model_name != NULL) {
        printf("model name provided, execution on model %s\n",
               options->model_name);

        instance inst = parse_input_file(options->model_name, "tsp", TSPLIB);
        add_params(inst, params);

        print_instance(inst, 1);
        solution sol = TSPopt(inst, BENDERS_CALLBACK);
        print_solution(sol, 1);
        plot_solution_graphviz(sol);

        return EXIT_SUCCESS;
    }

    /* execute battery test */
    options->battery_test = maxi(1, options->battery_test);
    printf("generating %d instances\n", options->battery_test);

    enum model_types tests[] = {
        MTZ_STATIC, MTZ_LAZY, GG_LAZY,
        BENDERS,  //   SYMMETRIC_BENDERS_CALLBACK,*/
    };

    int num_nodes = 20;
    instance* insts =
        generate_random_instances(options->battery_test, num_nodes, 20.0);

    for (int i = 0; i < options->battery_test; i++) {
        printf("battery %d:\n", i+1);

        instance inst = insts[i];
        add_params(inst, params);

        for (int j = 0; j < 4; j++) {
            char* model_type_str = model_type_tostring(tests[j]);
            printf("\tsolving %s on instance %s\n", model_type_str,
                   inst->model_name);
            free(model_type_str);

            TSPopt(inst, tests[j]);
        }
    }
    save_results(insts, options->battery_test);

    return EXIT_SUCCESS;
}
