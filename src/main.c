#include <assert.h>
#include <concorde.h>
#include <cplex.h>
#include <stdio.h>
#include <stdlib.h>

#include "../include/globals.h"
#include "../include/parsers.h"
#include "../include/solvers.h"
#include "../include/tsp.h"
#include "../include/utils.h"

int main(int argc, char** argv) {
    int ntests = 2;
    enum model_types tests[] = {/*MTZ_STATIC,*/
                                /*MTZ_LAZY,*/
                                /*GGLIT_STATIC,*/
                                /*GGLECT_STATIC,*/
                                /*GGLIT_LAZY,*/
                                /* BENDERS, */
                                BENDERS_CALLBACK,
                                HARD_FIXING};

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
            solution sol = TSPopt(inst, HARD_FIXING);
            print_solution(sol, 1);
            break;
        }
        case GENERATE: {
            options->battery_test = maxi(1, options->battery_test);
            printf("generating %d instances\n", options->battery_test);

            int num_nodes = 10;
            instance* insts = generate_random_instances(options->battery_test,
                                                        num_nodes, 20.0);
            double zstar = 0.0;

            for (int i = 0; i < options->battery_test; i++) {
                printf("battery %d:\n", i + 1);

                instance inst = insts[i];
                add_params(inst, params);

                for (int j = 0; j < ntests; j++) {
                    char* model_type_str = model_type_tostring(tests[j]);
                    printf("\tsolving %s on instance %s: ", model_type_str,
                           inst->model_name);
                    free(model_type_str);

                    solution sol = TSPopt(inst, tests[j]);
                    printf("%lf, ", sol->zstar);
                    printf("%lf\n", sol->solve_time);

                    if (zstar == 0.0) zstar = sol->zstar;
                    assert(fabs(zstar - sol->zstar) < EPSILON);
                }
            }
            save_results(insts, options->battery_test);

            break;
        }
        case LOAD_DIR: {
            printf("loading from directory\n");

            int ninstances;
            instance* insts =
                parse_input_dir(options->folder, "tsp", &ninstances, 0, 40);
            double zstar = 0.0;

            for (int i = 0; i < ninstances; i++) {
                instance inst = insts[i];
                printf("instance %s:\n", inst->model_name);

                add_params(inst, params);

                for (int j = 0; j < ntests; j++) {
                    char* model_type_str = model_type_tostring(tests[j]);
                    printf("\tsolving %s on instance %s: ", model_type_str,
                           inst->model_name);
                    free(model_type_str);

                    solution sol = TSPopt(inst, tests[j]);
                    printf("%lf, ", sol->zstar);
                    printf("%lf\n", sol->solve_time);

                    if (zstar == 0.0) zstar = sol->zstar;
                    assert(fabs(zstar - sol->zstar) < EPSILON);
                }
            }
            save_results(insts, ninstances);

            break;
        }
    }

    return EXIT_SUCCESS;
}
