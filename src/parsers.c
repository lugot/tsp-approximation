#include "parsers.h"
#include "tsp.h"
#include "utils.h"
#include "globals.h"
#include <getopt.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>

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
enum instance_types instance_type_enumerator (char* section_param);
enum weight_types weight_type_enumerator(char* section_param);

void print_usage() {
    printf("Usage: ./<name_executable> [options]\n");
    printf("Options:\n");
    printf("  --verbose                 Display more informations\n");
    printf("  --model_name <model_name> Model name present in data subfolder\n");
    printf("  --time_limit <integer>    Max timelimit for CPLEX computation\n");
    printf("  --seed <integer>          Randomness seed used in computation\n");
    printf("  --threads <integer>       Max number of threads\n");
    printf("  --memory <integer>        Max memory in MB used by CPLEX computation\n");
    printf("  --integer_costs           Consider integer costs only\n");
}

void parse_command_line(int argc, char** argv, instance inst) {

    static struct option long_options[] = {
        {"verbose",       no_argument,       NULL, 'v'},
        {"model_name",    required_argument, NULL, 'm'},
        {"time_limit",    required_argument, NULL, 't'},
        {"seed",          required_argument, NULL, 's'},
        {"threads",       required_argument, NULL, 'T'},
        {"memory",        required_argument, NULL, 'M'},
        {"integer_costs", no_argument,       NULL, 'i'},
        {"help",          no_argument,       NULL, 'h'},
        {0,               0,                 NULL, 0  }
    };

    int long_index, opt;
    long_index = opt = 0;
    while ((opt = getopt_long(argc, argv,"vm:t:s:T:M:i:h",
                    long_options, &long_index)) != -1) {
        switch (opt) {
            case 'v': VERBOSE = 1;
                break;

            case 'm':
                inst->model_name = (char*) malloc(strlen(optarg)*sizeof(char));
                strcpy(inst->model_name, optarg);
                break;

            case 'l': inst->timelimit = atof(optarg);
                break;

            case 's': inst->randomseed = abs(atoi(optarg));
                break;

            case 'T': inst->num_threads = atoi(optarg);
                break;

            case 'M': inst->available_memory = atoi(optarg);
                break;

            case 'c': inst->costs_type = INTEGER;
                break;

            case '?':
                print_usage();
                print_error("Unknown option `-%c'.\n", optopt);
        }
    }
}

void parse_input_file(instance inst, const char* file_extension) {

	if (VERBOSE) printf("[VERBOSE] reading input file %s\n", inst->model_name);

    char* filename;
    filename = (char*) calloc(100, sizeof(char));
    sprintf(filename, "../data/%s/%s.", inst->model_name, inst->model_name);
    strcat(filename, file_extension);

	FILE *fp;
    fp = fopen(filename, "r");
	if (fp == NULL) print_error("File %s not found", filename);


	char* line;
    char *section_name, *section_param;
    size_t len = 0;

	while((getline(&line, &len, fp) != -1)) {
		if (strlen(line) <= 1) continue;
        line[strlen(line)-1] = 0;  /* get rid of linefeed */

        section_name = strtok(line, " :");
        section_param = strtok(NULL, " :");

        if (VERBOSE) {
            printf("[VERBOSE] parsing |%s| on |%s|\n", section_name, section_param);
        }

        /* retrive section and inject parameter */
        switch (section_enumerator(section_name)) {

            /* specification part */
            case NAME:
                /*if we open the optimal tour file, this will overwrite its content,
                 * that makes incoherent the model name in common sense*/

                /*free(inst->model_name);*/
                /*inst->model_name = (char*) malloc(strlen(section_param)*sizeof(char));*/
                /*strcpy(inst->model_name, section_param);*/
                break;

            case TYPE:
                inst->instance_type = instance_type_enumerator(section_param);
                if (inst->instance_type == UNHANDLED_INSTANCE_TYPE) {
                    print_error("type %s unmanaged", section_param);
                }
                break;

            case COMMENT:
                inst->model_comment = (char*) malloc(strlen(section_param)*sizeof(char));
                strcpy(inst->model_comment, section_param);
                break;

            case DIMENSION:
                inst->num_nodes = atoi(section_param);
                break;

            case EDGE_WEIGHT_TYPE:
                inst->weight_type = weight_type_enumerator(section_param);
                if (inst->weight_type == UNHANDLED_WEIGHT_TYPE) {
                    print_error("weight type %s unhandled", section_param);
                }
                break;

            case EDGE_WEIGHT_FORMAT:
                if (strcmp(section_param, "LOWER_DIAG_ROW")) {
                    print_error("weight format %s unhandled", section_param);
                };
                break;


            /* data part */
            case NODE_COORD_SECTION:
            case DISPLAY_DATA_SECTION: {

                inst->nodes = (node*) calloc(inst->num_nodes, sizeof(node));

                int i;
                for (i=1; i<=inst->num_nodes && getline(&line, &len, fp) != -1; i++) {

                    int node_idx;
                    double x, y;

                    /*while (line[0] == ' ') line++;  [> avoid initial formatting spaces <]*/
                    node_idx = atoi(strtok(line, " "));
                    x = atof(strtok(NULL, " "));
                    y = atof(strtok(NULL, " "));

                    if(node_idx != i) print_error("incoherent node indexing\n");

                    inst->nodes[i-1] = (node) {x, y};
                }
                if (i-1 != inst->num_nodes) print_error("reached eof while reading nodes\n");

                break;
            }
            case TOUR_SECTION: {

                solution sol = (solution) calloc(1, sizeof(struct solution_t));
                sol->model_type = OPTIMAL_TOUR;
                sol->zstar = DBL_MAX;

                sol->num_edges = inst->num_nodes;
                sol->edges = (edge*) calloc(sol->num_edges, sizeof(struct edge_t));
                sol->parent = (int*) calloc(sol->num_edges, sizeof(int));

                int prev, first, act;
                if (getline(&line, &len, fp) == -1) print_error("reached eof while reading tour \n");
                prev = first = atoi(line)-1;

                int i;
                for (i=1; i<=inst->num_nodes-1 && getline(&line, &len, fp) != -1; i++) {
                    act = atoi(line)-1;

                    sol->parent[prev] = act;
                    sol->edges[i-1] = (edge) {prev, act};
                    prev = act;
                }
                sol->parent[prev] = first;
                sol->edges[i-1] = (edge) {prev, first};

                if (i != inst->num_nodes) print_error("reached eof while reading tour\n");

                add_solution(inst, sol);

                break;
            }

            case EDGE_WEIGHT_SECTION: {

                int num_weights = inst->num_nodes + inst->num_nodes*(inst->num_nodes-1)/2;
                double* full_weights = (double*) calloc(num_weights, sizeof(double));

                int k=0;
                while (getline(&line, &len, fp) && strcmp(line, "EOF\n")) {
                    char* weight_str;
                    weight_str = strtok(line, " ");

                    full_weights[k++] = atof(weight_str);

                    while ((weight_str = strtok(NULL, " ")) != NULL) {
                        full_weights[k++] = atof(weight_str);
                    }
                }

                k=0;
                inst->adjmatrix = (double**) calloc(inst->num_nodes, sizeof(double*));
                for (int i=0; i<inst->num_nodes; i++) {
                    inst->adjmatrix[i] = (double*) calloc(i+1, sizeof(double));

                    for (int j=0; j<=i; j++) inst->adjmatrix[i][j] = full_weights[k++];
                }

                break;
            }

            default:
            case END_OF_FILE:
                break;
            case UNHANDLED_SECTION:
                printf("Warning: section %s unmanaged\n", section_name);
                break;
        }
	}
    free(line);
    fclose(fp);
    free(filename);
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

    for (int i=NAME; i<=EDGE_WEIGHT_SECTION; i++) {
        if (!strcmp(section_name, sections[i])) return i;
    }

    return UNHANDLED_SECTION;
}
enum instance_types instance_type_enumerator (char* section_param) {
    char* instance_types[] = {
        "TSP",
        "TOUR",
    };

    for (int i=TSP; i<=TOUR; i++) {
        if (!strcmp(section_param, instance_types[i])) return i;
    }
    return UNHANDLED_INSTANCE_TYPE;
}
enum weight_types weight_type_enumerator(char* section_param) {
    char* weight_types[] = {
        "ATT",
        "EUC_2D",
        "GEO",
        "EXPLICIT"
    };

    for (int i=ATT; i<=EXPLICIT; i++) {
        if (!strcmp(section_param, weight_types[i])) return i;
    }
    return UNHANDLED_WEIGHT_TYPE;
}
