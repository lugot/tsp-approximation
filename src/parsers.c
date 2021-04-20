#include "../include/parsers.h"

#include <assert.h>
#include <float.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../include/globals.h"
#include "../include/tsp.h"
#include "../include/utils.h"

enum sections {
    NAME,
    TYPE,
    COMMENT,
    DIMENSION,
    CAPACITY,
    EDGE_WEIGHT_TYPE,
    EDGE_WEIGHT_FORMAT,
    EDGE_DATA_FORMAT,
    NODE_COORD_TYPE,
    DISPLAY_DATA_TYPE,
    END_OF_FILE,
    NODE_COORD_SECTION,
    DEPOT_SECTION,
    DEMAND_SECTION,
    EDGE_DATA_SECTION,
    FIXED_EDGE_SECTION,
    DISPLAY_DATA_SECTION,
    TOUR_SECTION,
    EDGE_WEIGHT_SECTION,
    UNHANDLED_SECTION,
};
enum sections section_enumerator(char* section_name);
enum instance_types instance_type_enumerator(char* section_param);
enum weight_types weight_type_enumerator(char* section_param);

run_options create_options() {
    run_options options = (run_options)calloc(1, sizeof(struct run_options_t));

    return options;
}

void print_usage() {
    printf("Usage: ./<name_executable> [options]\n");
    printf("Options:\n");
    printf("  --verbose                 Display more informations\n");
    printf(
        "  --model_name <model_name> Model name present in data subfolder\n");
    printf("  --time_limit <integer>    Max timelimit for CPLEX computation\n");
    printf("  --seed <integer>          Randomness seed used in computation\n");
    printf("  --threads <integer>       Max number of threads\n");
    printf(
        "  --memory <integer>        Max memory in MB used by CPLEX "
        "computation\n");
    printf("  --integer_costs           Consider integer costs only\n");
}

void parse_command_line(int argc, char** argv, cplex_params params,
                        run_options options) {
    assert(params != NULL);
    assert(options != NULL);

    static struct option long_options[] = {
        {"verbose", no_argument, NULL, 'v'},
        {"model_name", required_argument, NULL, 'm'},
        {"battery test", required_argument, NULL, 'b'},
        {"time_limit", required_argument, NULL, 't'},
        {"cplex_seed", required_argument, NULL, 's'},
        {"threads", required_argument, NULL, 'T'},
        {"memory", required_argument, NULL, 'M'},
        {"integer_costs", no_argument, NULL, 'i'},
        {"help", no_argument, NULL, 'h'},
        {0, 0, NULL, 0}};

    int long_index, opt;
    long_index = opt = 0;
    while ((opt = getopt_long(argc, argv, "vm:b:t:s:T:M:i:h", long_options,
                              &long_index)) != -1) {
        switch (opt) {
            case 'v':
                VERBOSE = 1;
                break;

            case 'm':
                options->model_name =
                    (char*)calloc(strlen(optarg), sizeof(char));
                snprintf(options->model_name, 1 + strlen(optarg), "%s", optarg);
                break;

            case 'b':
                options->battery_test = atoi(optarg);
                break;

            case 'l':
                params->timelimit = atof(optarg);
                break;

            case 's':
                params->randomseed = abs(atoi(optarg));
                break;

            case 'T':
                params->num_threads = atoi(optarg);
                break;

            case 'M':
                params->available_memory = atoi(optarg);
                break;

            case 'c':
                params->cost = INTEGER;
                break;

            case '?':
                print_usage();
                assert(opt != '?' && "unknown option");
        }
    }
}

instance parse_input_file(char* model_name, char* file_extension,
                          enum model_folders folder) {
    /* create instance from params */
    instance inst = create_empty_instance();

    assert(model_name != NULL &&
           "dont know what to parse, wrong execution mode");

    inst->model_name = (char*)calloc(strlen(model_name), sizeof(char));
    snprintf(inst->model_name, 1+strlen(model_name), "%s", model_name);

    inst->model_folder = folder;

    if (VERBOSE) printf("[VERBOSE] reading input file %s\n", inst->model_name);

    /* build filename depending on folder, model name and extention (tour or
     * tsp) */
    char* fname;
    fname = (char*)calloc(100, sizeof(char));
    if (folder == TSPLIB)
        snprintf(fname, strlen(fname), "../data_tsplib/%s/%s.",
                 inst->model_name, inst->model_name);
    if (folder == GENERATED)
        snprintf(fname, 17 + 2 * strlen(model_name), "../data_generated/%s/%s.",
                 inst->model_name, inst->model_name);
    snprintf(fname + strlen(fname), strlen(file_extension), "%s",
             file_extension);

    FILE* fp;
    fp = fopen(fname, "r");
    printf("filename: %s\n", fname);
    assert(fp != NULL && "file not found while parsing");

    char* line;
    char *section_name, *section_param;
    size_t len = 0;
    char *saveptr1, *saveptr2;

    while ((getline(&line, &len, fp) != -1)) {
        if (strlen(line) <= 1) continue;
        line[strlen(line) - 1] = 0; /* get rid of linefeed */

        section_name = strtok_r(line, " :", &saveptr1);
        section_param = strtok_r(NULL, " :", &saveptr2);

        if (VERBOSE)
            printf("[Verbose] parsing |%s| on |%s|\n", section_name,
                   section_param);

        /* retrive section and inject parameter */
        switch (section_enumerator(section_name)) {
            /* specification part */
            case NAME:
                /*if we open the optimal tour file, this will overwrite its
                 * content, that makes incoherent the model name in common
                 * sense*/

                /*free(inst->model_name);*/
                /*inst->model_name = (char*)
                 * malloc(strlen(section_param)*sizeof(char));*/
                /*strcpy(inst->model_name, section_param);*/
                break;

            case TYPE:
                inst->instance_type = instance_type_enumerator(section_param);
                assert(inst->instance_type != UNHANDLED_INSTANCE_TYPE &&
                       "able to handly TPS or TOUR");
                break;

            case COMMENT:
                inst->model_comment =
                    (char*)malloc(strlen(section_param) * sizeof(char));
                snprintf(inst->model_comment, strlen(section_param), "%s",
                         section_param);
                break;

            case DIMENSION:
                inst->nnodes = atoi(section_param);
                break;

            case EDGE_WEIGHT_TYPE:
                inst->weight_type = weight_type_enumerator(section_param);
                assert(inst->weight_type != UNHANDLED_WEIGHT_TYPE &&
                       "weight type unhandled");
                break;

            case EDGE_WEIGHT_FORMAT:
                assert(strcmp(section_param, "LOWER_DIAG_ROW") == 0 &&
                       "edge weight format unhandled");
                break;

            /* data part */
            case NODE_COORD_SECTION:
            case DISPLAY_DATA_SECTION: {
                int nnodes = inst->nnodes;
                inst->nodes = (node*)calloc(nnodes, sizeof(struct node_t));

                int i;
                for (i = 1; i <= nnodes; i++) {
                    /* reached eof */
                    if (getline(&line, &len, fp) != -1) break;

                    int node_idx;
                    double x, y;

                    node_idx = atoi(strtok_r(line, " ", &saveptr2));
                    x = atof(strtok_r(NULL, " ", &saveptr2));
                    y = atof(strtok_r(NULL, " ", &saveptr2));

                    assert(node_idx == i &&
                           "incoherent node indexing in NODE_COORD or "
                           "DISPLAY_DATA sections");

                    inst->nodes[i - 1] = (node){x, y};
                }
                assert(i - 1 == nnodes && "reached EOF while reading nodes");

                break;
            }
            case TOUR_SECTION: {
                /* a tour is managed by creating an empty instance with a single
                 * solution (the tour) */
                solution sol = (solution)calloc(1, sizeof(struct solution_t));
                sol->model_type = OPTIMAL_TOUR;
                /* the tour does not provide the min cost */
                sol->zstar = DBL_MAX;

                int nedges = sol->nedges = inst->nnodes;
                sol->edges = (edge*)calloc(nedges, sizeof(struct edge_t));
                sol->link = (int*)calloc(nedges, sizeof(int));

                int prev, first;
                assert(getline(&line, &len, fp) != -1 &&
                       "reached EOF of while reading tour");
                prev = first = atoi(line) - 1;

                int i;
                for (i = 1; i <= inst->nnodes - 1; i++) {
                    if (getline(&line, &len, fp) != -1) break;
                    int act = atoi(line) - 1;

                    sol->link[prev] = act;
                    sol->edges[i - 1] = (edge){prev, act};
                    prev = act;
                }
                sol->link[prev] = first;
                sol->edges[i - 1] = (edge){prev, first};

                assert(i == inst->nnodes && "reached EOF while reading tour");

                add_solution(inst, sol);

                break;
            }

            case EDGE_WEIGHT_SECTION: {
                int nnodes = inst->nnodes;
                int nweights = nnodes + nnodes * (nnodes - 1) / 2;
                double* weights = (double*)calloc(nweights, sizeof(double));

                int k = 0;
                while (getline(&line, &len, fp) && strcmp(line, "EOF\n")) {
                    char* weight_str;
                    weight_str = strtok_r(line, " ", &saveptr2);

                    weights[k++] = atof(weight_str);

                    while ((weight_str = strtok_r(NULL, " ", &saveptr2)) !=
                           NULL) {
                        weights[k++] = atof(weight_str);
                    }
                }

                k = 0;
                inst->adjmatrix = (double**)calloc(nnodes, sizeof(double*));
                for (int i = 0; i < nnodes; i++) {
                    inst->adjmatrix[i] = (double*)calloc(i + 1, sizeof(double));

                    for (int j = 0; j <= i; j++)
                        inst->adjmatrix[i][j] = weights[k++];
                }

                break;
            }

            default:
            case END_OF_FILE:
                break;
            case UNHANDLED_SECTION:
                if (VERBOSE)
                    printf("[Warning] section %s unmanaged\n", section_name);
                break;
        }
    }
    free(line);
    fclose(fp);
    free(fname);

    return inst;
}

enum sections section_enumerator(char* section_name) {
    char* sections[] = {
        "NAME",
        "TYPE",
        "COMMENT",
        "DIMENSION",
        "CAPACITY",
        "EDGE_WEIGHT_TYPE",
        "EDGE_WEIGHT_FORMAT",
        "EDGE_DATA_FORMAT",
        "NODE_COORD_TYPE",
        "DISPLAY_DATA_TYPE",
        "EOF",
        "NODE_COORD_SECTION",
        "DEPOT_SECTION",
        "DEMAND_SECTION",
        "EDGE_DATA_SECTION",
        "FIXED_EDGE_SECTION",
        "DISPLAY_DATA_SECTION",
        "TOUR_SECTION",
        "EDGE_WEIGHT_SECTION",
    };

    for (int i = NAME; i <= EDGE_WEIGHT_SECTION; i++) {
        if (!strcmp(section_name, sections[i])) return i;
    }

    return UNHANDLED_SECTION;
}
enum instance_types instance_type_enumerator(char* section_param) {
    char* instance_types[] = {
        "TSP",
        "TOUR",
    };

    for (int i = TSP; i <= TOUR; i++) {
        if (!strcmp(section_param, instance_types[i])) return i;
    }
    return UNHANDLED_INSTANCE_TYPE;
}
enum weight_types weight_type_enumerator(char* section_param) {
    char* weight_types[] = {"ATT", "EUC_2D", "GEO", "EXPLICIT"};

    for (int i = ATT; i <= EXPLICIT; i++) {
        if (!strcmp(section_param, weight_types[i])) return i;
    }
    return UNHANDLED_WEIGHT_TYPE;
}
