#include "tsp.h"
#include "parsers.h"
#include "utils.h"
#include "solvers.h"
#include <stdio.h>
#include <cplex.h>
#include <stdlib.h>
#include <assert.h>

int main(int argc, char** argv) {

    cplex_params params = (cplex_params) calloc(1, sizeof(struct cplex_params_t));
    run_options options = (run_options) calloc(1, sizeof(struct run_options_t));
    parse_command_line(argc, argv, params, options);

    if (options->model_name != NULL) {
        printf("model name provided, execution on model %s\n", options->model_name);

        instance inst = parse_input_file(options->model_name, "tsp", TSPLIB);
        add_params(inst, params);

        print_instance(inst, 1);
        solution sol = TSPopt(inst, SYMMETRIC_BENDERS_CALLBACK);
        print_solution(sol, 1);
        plot_solution_graphviz(sol);

        return EXIT_SUCCESS;
    }

    /* execute battery test */
    options->battery_test = maxi(1, options->battery_test);
    printf("generating %d instances\n", options->battery_test);

    enum model_types tests[] = {
        SYMMETRIC_BENDERS,
        SYMMETRIC_BENDERS_CALLBACK
    };

    int num_nodes = 20;
    instance* insts = generate_random_instances(options->battery_test, num_nodes, 20.0);

    for (int i=0; i<options->battery_test; i++) {
        instance inst = insts[i];
        add_params(inst, params);

        for (int i=0; i<2; i++) {
            char* model_type_str = model_type_tostring(tests[i]);
            printf("solving %s on instance %s\n", model_type_str, inst->model_name);
            free(model_type_str);

            TSPopt(inst, tests[i]);
        }
        /*print_solution(sol, 1);*/
        /*plot_solution_graphviz(sol);*/
        /*print_instance(inst, 1);*/
    }
    save_results(insts, options->battery_test);


    return EXIT_SUCCESS;
}
