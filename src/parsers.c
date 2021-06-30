#include "../include/parsers.h"

#include <assert.h>
#include <errno.h>
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
    options->mode = NOT_SPECIFIED;

    return options;
}
void free_options(run_options options) {
    free(options->model_name);
    free(options);
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
    printf("  models:\n");
    for (int i = 0; i < 23; i++) {
        char* model_type_str = model_type_tostring(i);
        printf("\t%s: %d %d\n", model_type_str, 1 << i, (1 << (i + 1)) - 1);
        free(model_type_str);
    }
}

void parse_command_line(int argc, char** argv, cplex_params params,
                        run_options options) {
    assert(params != NULL);
    assert(options != NULL);

    static struct option long_options[] = {
        {"verbose", no_argument, NULL, 'V'},
        {"extra", no_argument, NULL, 'E'},
        {"suppress", no_argument, NULL, 'S'},
        {"gen_nnodes", no_argument, NULL, 'N'},
        {"model_name", required_argument, NULL, 'm'},
        {"generate", required_argument, NULL, 'g'},
        {"load_directory", required_argument, NULL, 'l'},
        {"tests", required_argument, NULL, 'W'},
        {"optimal", required_argument, NULL, 'O'},
        {"time_limit", required_argument, NULL, 't'},
        {"cplex_seed", required_argument, NULL, 's'},
        {"threads", required_argument, NULL, 'T'},
        {"memory", required_argument, NULL, 'M'},
        {"integer_costs", no_argument, NULL, 'i'},
        {"help", no_argument, NULL, 'h'},
        {0, 0, NULL, 0}};

    int long_index, opt;
    long_index = opt = 0;
    while ((opt = getopt_long(argc, argv, "VESON:m:g:l:W:b:t:s:T:M:i:h",
                              long_options, &long_index)) != -1) {
        switch (opt) {
            case 'V':
                VERBOSE = 1;
                break;
            case 'E':
                EXTRA = 1;
                break;
            case 'S':
                SUPPRESS_CALLBACK = 0;
                break;
            case 'N':
                GEN_NNODES = atoi(optarg);
                assert(GEN_NNODES > 5);
                break;
            case 'm':
                options->mode = SINGLE_MODEL;
                options->model_name =
                    (char*)calloc(1 + strlen(optarg), sizeof(char));
                snprintf(options->model_name, 1 + strlen(optarg), "%s", optarg);
                break;
            case 'g':
                options->mode = GENERATE;
                options->battery_test = atoi(optarg);
                break;
            case 'l':
                options->mode = LOAD_DIR;
                if (!strcmp(optarg, "tsplib"))
                    options->folder = TSPLIB;
                else if (!strcmp(optarg, "generated"))
                    options->folder = GENERATED;
                else
                    print_error("wrong folder");
                break;
            case 'W':
                options->tests = atoi(optarg);
                assert(GEN_NNODES > 5);
                break;
            case 'O':
                options->load_optimal = 1;
                break;
            case 't':
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
                printf("opt: %c, optarg: %s\n", opt, optarg);
                assert(opt != '?' && "unknown option");
        }
    }
    if (options->mode == NOT_SPECIFIED) {
        print_usage();
        printf("opt: %c, optarg: %s\n", opt, optarg);
        assert(options->mode != NOT_SPECIFIED && "use a mode");
    }
}

instance parse_input_file(char* model_name, char* file_extension,
                          enum model_folders folder) {
    /* create instance from params */
    instance inst = create_empty_instance();

    assert(model_name != NULL &&
           "dont know what to parse, wrong execution mode");

    inst->model_name = (char*)calloc(1 + strlen(model_name), sizeof(char));
    snprintf(inst->model_name, 1 + strlen(model_name), "%s", model_name);

    inst->model_folder = folder;

    if (VERBOSE) printf("[VERBOSE] reading input file %s\n", inst->model_name);

    /* build filename depending on folder, model name and extention (tour or
     * tsp) */
    char* fname;
    int bufsize = 100;
    fname = (char*)calloc(bufsize, sizeof(char));
    if (folder == TSPLIB)
        snprintf(fname, bufsize, "../data_tsplib/%s/%s.", inst->model_name,
                 inst->model_name);
    if (folder == GENERATED)
        snprintf(fname, bufsize, "../data_generated/%s/%s.", inst->model_name,
                 inst->model_name);
    snprintf(fname + strlen(fname), bufsize - strlen(fname), "%s",
             file_extension);

    FILE* fp;
    fp = fopen(fname, "r");
    assert(fp != NULL && "file not found while parsing");

    char* line = "";
    char *section_name, *section_param;
    size_t len = 0;
    char* saveptr = "";

    while ((getline(&line, &len, fp) != -1)) {
        if (strlen(line) <= 1) continue;
        line[strlen(line) - 1] = 0; /* get rid of linefeed */

        int icolon = 0;
        while (icolon < strlen(line) && line[icolon] != ':') icolon++;
        if (icolon == strlen(line)) {
            /* parsing NODE_COORD_SECTION or similar, no param */
            section_name = line;
            section_param = NULL;
        } else if (line[icolon - 1] == ' ') {
            /* parsing att48-like tsp files */
            section_name = strtok_r(line, " : ", &saveptr);
            section_param = strtok_r(NULL, "\n", &saveptr);
            section_param += 2;
        } else {
            /* parsing berlin52-like tsp files */
            section_name = strtok_r(line, ": ", &saveptr);
            section_param = strtok_r(NULL, "\n", &saveptr);
            section_param += 1;
        }

        if (VERBOSE) {
            printf("[Verbose] parsing %s, %s\n", section_name, section_param);
        }

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
                    (char*)calloc(1 + strlen(section_param), sizeof(char));
                snprintf(inst->model_comment, 1 + strlen(section_param), "%s",
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
                assert(strcmp(section_param, "LOWER_DIAG_ROW") != 0 &&
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
                    if (getline(&line, &len, fp) == -1) break;

                    int node_idx;
                    double x, y;

                    node_idx = atoi(strtok_r(line, " ", &saveptr));
                    x = atof(strtok_r(NULL, " ", &saveptr));
                    y = atof(strtok_r(NULL, " ", &saveptr));

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

                int prev, first;
                assert(getline(&line, &len, fp) != -1 &&
                       "reached EOF of while reading tour");
                prev = first = atoi(line) - 1;

                int i;
                for (i = 1; i <= inst->nnodes - 1; i++) {
                    if (getline(&line, &len, fp) == -1) break;
                    int act = atoi(line) - 1;

                    sol->edges[i - 1] = (edge){prev, act};
                    prev = act;
                }
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
                    weight_str = strtok_r(line, " ", &saveptr);

                    weights[k++] = atof(weight_str);

                    while ((weight_str = strtok_r(NULL, " ", &saveptr)) !=
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
                free(weights);

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

instance* parse_input_dir(enum model_folders folder, char* file_extension,
                          int* ninstances, int nodes_lb, int nodes_ub) {
    assert(nodes_lb >= 0);
    assert(nodes_ub > 0);
    assert(nodes_lb < nodes_ub);
    char** model_names = list_files(folder, ninstances);

    instance* insts = (instance*)calloc(*ninstances, sizeof(struct instance_t));

    int nmodels = *ninstances;
    *ninstances = 0;
    for (int i = 0; i < nmodels; i++) {
        instance inst =
            parse_input_file(model_names[i], file_extension, folder);

        if (nodes_lb <= inst->nnodes && inst->nnodes < nodes_ub) {
            insts[*ninstances] = inst;
            *ninstances = *ninstances + 1;
        }

        free(model_names[i]);
    }
    free(model_names);

    instance* tmp =
        (instance*)realloc(insts, (*ninstances) * sizeof(struct instance_t));
    if (tmp == NULL) {
        free(insts);
        assert(tmp != NULL);

    } else {
        insts = tmp;
    }

    return insts;
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
