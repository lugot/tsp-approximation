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
    options->tests = -1;

    return options;
}
void free_options(run_options options) {
    free(options->instance_name);
    free(options->instance_folder);

    free(options);
}

void print_usage() {
    printf("Usage: ./<name_executable> [options]\n");
    printf("  -v --verbose\n");
    printf("  -e --extra_verbose\n");
    printf("  -c --callback_verbose\n");
    printf("  -n --instance_name <instance_name>\n");
    printf("  -l --load_directory <subfolder in ./data>\n");
    printf("  -o --load_optimal\n");
    printf("  -g --generate <number of instance to generate>\n");
    printf("  -N --nnodes <number of nodes of generated instances>\n");
    printf("  -m --models <execute models>\n");
    printf("  -T --time_limit <time limit in seconds>\n");
    printf("  -S --cplex_seed <cplex seed>\n");
    printf("  -C --threads <threads to use>\n");
    printf("  -M --memory <max memory usage in MB>\n");
    printf("  -h --help\n");
    printf("  avaiable models:\n");
    for (int i = 0; i < 24; i++) {
        char* model_type_str = model_type_tostring(i);
        printf("\t%s: %d\n", model_type_str, 1 << i);
        free(model_type_str);
    }

    exit(EXIT_SUCCESS);
}

void parse_command_line(int argc, char** argv, cplex_params params,
                        run_options options) {
    assert(params != NULL);
    assert(options != NULL);

    static struct option long_options[] = {
        {"verbose", no_argument, NULL, 'v'},
        {"extra_verbose", no_argument, NULL, 'e'},
        {"callback_verbose", no_argument, NULL, 'c'},
        {"instance_name", required_argument, NULL, 'n'},
        {"load_directory", required_argument, NULL, 'l'},
        {"load_optimal", no_argument, NULL, 'o'},
        {"generate", required_argument, NULL, 'g'},
        {"nnodes", required_argument, NULL, 'N'},
        {"models", required_argument, NULL, 'm'},
        {"time_limit", required_argument, NULL, 'T'},
        {"cplex_seed", required_argument, NULL, 'S'},
        {"threads", required_argument, NULL, 'C'},
        {"memory", required_argument, NULL, 'M'},
        {"help", no_argument, NULL, 'h'},
        {0, 0, NULL, 0}};

    int long_index, opt;
    long_index = opt = 0;
    while ((opt = getopt_long(argc, argv, "vecn:l:og:N:m:T:S:C:M:h",
                              long_options, &long_index)) != -1) {
        switch (opt) {
            case 'v':
                VERBOSE = 1;
                break;
            case 'e':
                EXTRA_VERBOSE = 1;
                break;
            case 'c':
                CALLBACK_VERBOSE = 1;
                break;
            case 'n':
                if (options->mode != NOT_SPECIFIED) {
                    print_error("cannot mix execution modes!");
                }
                options->mode = SINGLE_INSTANCE;
                options->instance_name =
                    (char*)calloc(1 + strlen(optarg), sizeof(char));
                snprintf(options->instance_name, 1 + strlen(optarg), "%s",
                         optarg);
                break;
            case 'l':
                if (options->mode != NOT_SPECIFIED) {
                    print_error("cannot mix execution modes!");
                }
                options->mode = LOAD_DIR;
                options->instance_folder =
                    (char*)calloc(1 + strlen(optarg), sizeof(char));
                snprintf(options->instance_folder, 1 + strlen(optarg), "%s",
                         optarg);
                break;
            case 'o':
                options->load_optimal = 1;
                break;
            case 'g':
                if (options->mode != NOT_SPECIFIED) {
                    print_error("cannot mix execution modes!");
                }
                options->mode = GENERATE;
                options->battery_test = atoi(optarg);
                break;
            case 'N':
                GEN_NNODES = atoi(optarg);
                assert(GEN_NNODES > 3);
                break;
            case 'm':
                options->tests = atoi(optarg);
                break;
            case 'T':
                params->timelimit = atof(optarg);
                break;
            case 'S':
                params->randomseed = abs(atoi(optarg));
                break;
            case 'C':
                params->num_threads = atoi(optarg);
                break;
            case 'M':
                params->available_memory = atoi(optarg);
                break;
            case 'h':
                print_usage();
                break;
            case '?':
                printf("opt: %c, optarg: %s\n", opt, optarg);
                print_usage();
        }
    }
    if (options->mode == NOT_SPECIFIED ||
        (options->mode != GENERATE && options->tests == -1)) {
        print_usage();
    }
}

instance parse_input_file(char* instance_name, char* file_extension,
                          char* instance_folder) {
    /* create instance from params */
    instance inst = create_empty_instance();

    assert(instance_name != NULL &&
           "dont know what to parse, wrong execution mode");

    inst->instance_name =
        (char*)calloc(1 + strlen(instance_name), sizeof(char));
    snprintf(inst->instance_name, 1 + strlen(instance_name), "%s",
             instance_name);

    inst->instance_folder =
        (char*)calloc(1 + strlen(instance_folder), sizeof(char));
    snprintf(inst->instance_folder, 1 + strlen(instance_folder), "%s",
             instance_folder);

    /* build filename depending on folder, model name and extention (tour or
     * tsp) */
    char* fname;
    int bufsize = 100;
    fname = (char*)calloc(bufsize, sizeof(char));

    snprintf(fname, bufsize, "../data/%s/%s/%s.%s", instance_folder,
             instance_name, instance_name, file_extension);

    if (VERBOSE) {
        printf("[VERBOSE] reading input file %s\n", instance_name);
    }

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
                /* not handled */
                break;

            case COMMENT:
                inst->instance_comment =
                    (char*)calloc(1 + strlen(section_param), sizeof(char));
                snprintf(inst->instance_comment, 1 + strlen(section_param),
                         "%s", section_param);
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
                solution sol =
                    create_solution(inst, OPTIMAL_TOUR, inst->nnodes);
                /* the tour does not provide the min cost */
                sol->zstar = DBL_MAX;

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

instance* parse_input_dir(char* instance_folder, char* file_extension,
                          int* ninstances, int load_optimal) {
    return parse_input_dir_bounded(instance_folder, file_extension, ninstances,
                                   load_optimal, 0, 10000000);
}
instance* parse_input_dir_bounded(char* instance_folder, char* file_extension,
                                  int* ninstances, int load_optimal,
                                  int nodes_lb, int nodes_ub) {
    assert(nodes_lb >= 0);
    assert(nodes_ub > 0);
    assert(nodes_lb < nodes_ub);
    char** instance_names = list_files(instance_folder, ninstances);

    instance* insts = (instance*)calloc(*ninstances, sizeof(struct instance_t));

    int nmodels = *ninstances;
    *ninstances = 0;
    for (int i = 0; i < nmodels; i++) {
        instance inst = parse_input_file(instance_names[i], file_extension,
                                         instance_folder);

        if (load_optimal) {
            instance opt = parse_input_file(instance_names[i], "opt.tour",
                                            instance_folder);

            inst->zbest = compute_zstar(inst, opt->sols[0]);
            free_instance(opt);
        }

        if (nodes_lb <= inst->nnodes && inst->nnodes < nodes_ub) {
            insts[*ninstances] = inst;
            *ninstances = *ninstances + 1;
        }

        free(instance_names[i]);
    }
    free(instance_names);

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
    };

    for (int i = NAME; i <= TOUR_SECTION; i++) {
        if (!strcmp(section_name, sections[i])) return i;
    }

    return UNHANDLED_SECTION;
}
enum weight_types weight_type_enumerator(char* section_param) {
    char* weight_types[] = {"ATT", "EUC_2D", "GEO", "EXPLICIT"};

    for (int i = ATT; i <= EXPLICIT; i++) {
        if (!strcmp(section_param, weight_types[i])) return i;
    }
    return UNHANDLED_WEIGHT_TYPE;
}
