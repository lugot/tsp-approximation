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
    cplex_params params = create_params();
    run_options options = create_options();
    parse_command_line(argc, argv, params, options);

    int ntests = 0;
    enum model_types tests[100];

    int testsint = options->tests;
    int i = 0;
    while (testsint > 0) {
        if (testsint % 2 == 1) tests[ntests++] = (enum model_types)i;

        i++;
        testsint /= 2;
    }

    printf("execution on:\n");
    for (int i = 0; i < ntests; i++) {
        char* model_type_str = model_type_tostring(tests[i]);
        printf("\t%s\n", model_type_str);
        free(model_type_str);
    }

    switch (options->mode) {
        case SINGLE_INSTANCE: {
            printf("model name provided, execution on model %s\n",
                   options->instance_name);

            instance inst =
                parse_input_file(options->instance_name, "tsp", "tsplib");
            add_params(inst, params);

            if (options->load_optimal) {
                instance opt = parse_input_file(inst->instance_name, "opt.tour",
                                                inst->instance_folder);

                inst->zbest = compute_zstar(inst, opt->sols[0]);
                free(opt);
            }

            print_instance(inst, 1);
            solution sol = solve(inst, tests[0]);

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
                generate_random_instances(options->battery_test, GEN_NNODES);

            for (int i = 0; i < options->battery_test && ntests != 0; i++) {
                add_params(insts[i], params);
                save_instance(insts[i]);

                if (EXTRA_VERBOSE) print_instance(insts[i], 1);
            }

            for (int i = 0; i < options->battery_test; i++) {
                free_instance(insts[i]);
            }
            free(insts);
            break;
        }
        case LOAD_DIR: {
            printf("loading from directory %s\n", options->instance_folder);

            FILE* emergency_res;
            emergency_res = fopen("emergency_res.csv", "a");
            assert(emergency_res != NULL && "file not found while saving emer");

            fprintf(emergency_res, "%d,", ntests);
            for (int j = 0; j < ntests; j++) {
                if (j == ntests - 1) {
                    char* model_type_str = model_type_tostring(tests[j]);
                    fprintf(emergency_res, "%s\n", model_type_str);
                    free(model_type_str);
                } else {
                    char* model_type_str = model_type_tostring(tests[j]);
                    fprintf(emergency_res, "%s,", model_type_str);
                    free(model_type_str);
                }
            }
            fclose(emergency_res);

            int ninstances;
            instance* insts =
                parse_input_dir(options->instance_folder, "tsp", &ninstances,
                                options->load_optimal);

            for (int i = 0; i < ninstances; i++) {
                instance inst = insts[i];
                add_params(inst, params);
                printf("instance %s:\n", inst->instance_name);

                FILE* emergency_res;
                emergency_res = fopen("emergency_res.csv", "a");
                assert(emergency_res != NULL &&
                       "file not found while saving emer");

                fprintf(emergency_res, "%s,", inst->instance_name);
                fclose(emergency_res);

                for (int j = 0; j < ntests; j++) {
                    /* solve! */
                    solution sol = solve(inst, tests[j]);

                    char* model_type_str = model_type_tostring(tests[j]);
                    printf("\tsolving %s on instance %s: %lf %lf\n",
                           model_type_str, inst->instance_name, sol->zstar,
                           sol->solve_time);
                    free(model_type_str);

                    FILE* emergency_res;
                    emergency_res = fopen("emergency_res.csv", "a");
                    assert(emergency_res != NULL &&
                           "file not found while saving emer");

                    if (j == ntests - 1) {
                        fprintf(emergency_res, "%lf\n", sol->solve_time);
                    } else {
                        fprintf(emergency_res, "%lf,", sol->solve_time);
                    }
                    fclose(emergency_res);

                    if (EXTRA_VERBOSE) plot_graphviz(sol, NULL, j);
                }
            }
            plot_profiler(insts, ninstances);
            if (options->load_optimal) plot_tracking(insts, ninstances, 50);
            if (options->load_optimal) plot_tracking(insts, ninstances, 10);
            if (options->load_optimal) plot_tracking(insts, ninstances, 3);

            for (int i = 0; i < ninstances; i++) {
                free_instance(insts[i]);
            }
            free(insts);

            break;
        }
        case NOT_SPECIFIED:
            break;
    }

    free_options(options);
    free_params(params);

    return EXIT_SUCCESS;
}
