#include <assert.h>
#include <concorde.h>
#include <cplex.h>
#include <stdio.h>
#include <stdlib.h>

#include "../include/globals.h"
#include "../include/parsers.h"
#include "../include/pqueue.h"
#include "../include/solvers.h"
#include "../include/tsp.h"
#include "../include/union_find.h"
#include "../include/utils.h"

int main(int argc, char** argv) {
    int ntests = 1;
    enum model_types tests[] = {/*MTZ_STATIC,*/
                                /*MTZ_LAZY,*/
                                /*GGLIT_STATIC,*/
                                /*GGLECT_STATIC,*/
                                /*GGLIT_LAZY,*/
                                /* BENDERS, */
                                /* BENDERS_CALLBACK, */
                                /* HARD_FIXING, */
                                /* SOFT_FIXING, */
                                /* MST, */
                                GREEDY,
                                /* GRASP, */
                                /* EXTRA_MILEAGE */
    };

    cplex_params params = create_params();
    run_options options = create_options();
    parse_command_line(argc, argv, params, options);

    switch (options->mode) {
        case SINGLE_MODEL: {
            printf("model name provided, execution on model %s\n",
                   options->model_name);

            instance inst =
                parse_input_file(options->model_name, "tsp", TSPLIB);
            add_params(inst, params);

            print_instance(inst, 1);
            solution sol = solve(inst, EXTRA_MILEAGE);

            /* single istance needs to be verbose */
            print_solution(sol, 1);
            plot_graphviz(sol, NULL, 0);

            free_instance(inst);
            break;
        }
        case GENERATE: {
            options->battery_test = maxi(1, options->battery_test);
            printf("generating %d instances\n", options->battery_test);

            /* generate and solve the instances */
            instance* insts =
                generate_random_instances(options->battery_test, GEN_NUM_NODES);
            for (int i = 0; i < options->battery_test; i++) {
                printf("battery %d:\n", i + 1);

                instance inst = insts[i];
                add_params(inst, params);
                save_instance(inst);

                for (int j = 0; j < ntests; j++) {
                    /* solve! */
                    solution sol = solve(inst, tests[j]);

                    char* model_type_str = model_type_tostring(tests[j]);
                    printf("\tsolving %s on instance %s: %lf %lf\n", model_type_str,
                           inst->model_name, sol->zstar, sol->solve_time);
                    free(model_type_str);

                    if (EXTRA) plot_graphviz(sol, NULL, j);
                }

                if (VERBOSE) print_instance(inst, 1);
            }
            plot_profiler(insts, options->battery_test);

            for (int i = 0; i < options->battery_test; i++) {
                free_instance(insts[i]);
            }
            free(insts);
            break;
        }
        case LOAD_DIR: {
            printf("loading from directory\n");

            int ninstances;
            instance* insts =
                parse_input_dir(options->folder, "tsp", &ninstances, 0, 40);

            for (int i = 0; i < ninstances; i++) {
                instance inst = insts[i];
                add_params(inst, params);
                printf("instance %s:\n", inst->model_name);

                for (int j = 0; j < ntests; j++) {
                    /* solve! */
                    solution sol = solve(inst, tests[j]);

                    char* model_type_str = model_type_tostring(tests[j]);
                    printf("\tsolving %s on instance %s: %lf %lf\n", model_type_str,
                           inst->model_name, sol->zstar, sol->solve_time);
                    free(model_type_str);

                    if (EXTRA) plot_graphviz(sol, NULL, j);
                }
            }
            plot_profiler(insts, ninstances);

            for (int i = 0; i < ninstances; i++) {
                free_instance(insts[i]);
            }
            free(insts);

            break;
        }
    }

    free_options(options);
    free_params(params);

    return EXIT_SUCCESS;
}
