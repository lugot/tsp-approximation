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

            for (int i = 0; i < options->battery_test; i++) {
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

            int ninstances;
            instance* insts =
                parse_input_dir(options->instance_folder, "tsp", &ninstances,
                                options->load_optimal);

            for (int i = 0; i < ninstances; i++) {
                instance inst = insts[i];
                add_params(inst, params);
                printf("instance %s:\n", inst->instance_name);

                for (int j = 0; j < ntests; j++) {
                    /* solve! */
                    solution sol = solve(inst, tests[j]);

                    char* model_type_str = model_type_tostring(tests[j]);
                    printf("\tsolving %s on instance %s: %lf %lf\n",
                           model_type_str, inst->instance_name, sol->zstar,
                           sol->solve_time);
                    free(model_type_str);

                    if (EXTRA_VERBOSE) plot_graphviz(sol, NULL, j);
                }

                /* plot_profiler(insts, i + 1, 0); */
                plot_profiler(insts, i + 1, 1);
                if (options->load_optimal) plot_tracking(insts, i + 1, 50);
                if (options->load_optimal) plot_tracking(insts, i + 1, 10);
                if (options->load_optimal) plot_tracking(insts, i + 1, 3);
            }

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
