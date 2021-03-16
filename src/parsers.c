#include "parsers.h"
#include "utils.h"
#include "globals.h"
#include <getopt.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

enum tsp_sections {
    NAME,
    COMMENT,
    TYPE,
    DIMENSION,
    EDGE_WEIGHT_TYPE,
    NODE_COORD_SECTION,
    TOUR_SECTION,
    END_OF_FILE,
    UNMANAGED_SECTION,
};
enum tsp_sections hash2section(char* section_name);

// TODO: add this to a help
void print_usage() {
    printf("Usage: rectangle -file num -b num\n");
}

// TODO: comment
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

            case 't': inst->model_type = SYMMETRIC_TSP;
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
                printf("Unknown option `-%c'.\n", optopt);
                print_usage();
                exit(EXIT_FAILURE);
        }
    }
}


// TODO: comment function
// TODO: change tokenizer
void parse_input_file(instance inst) {

	if (VERBOSE) printf("[VERBOSE] reading input file %s\n", inst->model_name);

    char* filename;
    filename = (char*) calloc(100, sizeof(char));
    sprintf(filename, "../data/%s/%s.tsp", inst->model_name, inst->model_name);

	FILE *fp;
    fp = fopen(filename, "r");
	if (fp == NULL) {
        char* errormessage;
        errormessage = (char*) calloc(100, sizeof(char));
        sprintf(errormessage, "File %s not found. Aborted\n", filename);
        print_error(errormessage);
    }


	char* line;
    char *section_name, *section_param;
    size_t len = 0;

	while((getline(&line, &len, fp) != -1)) {
		if (strlen(line) <= 1) continue;
        line[strlen(line)-1] = 0;  /* get rid of linefeed */

        /* old style tokenization */
        size_t i=0;
        while (line[i] != ' ') i++;

        section_name = (char*) calloc(i, sizeof(char));
        strncpy(section_name, line, i);
        i = (i+3 > strlen(line) ? strlen(line) : i+3);
        section_param = (char*) calloc(strlen(line)-i, sizeof(char));
        strcpy(section_param, line+i);

        if (VERBOSE) {
            printf("[VERBOSE] parsing |%s| on |%s|\n", section_name, section_param);
        }

        /* retrive section and inject parameter */
        switch (hash2section(section_name)) {
            case NAME:
                inst->model_name = (char*) malloc(strlen(section_param)*sizeof(char));
                strcpy(inst->model_name, section_param);
                break;

            case COMMENT:
                inst->model_comment = (char*) malloc(strlen(section_param)*sizeof(char));
                strcpy(inst->model_comment, section_param);
                break;

            case TYPE:
                /* keep the one set by user */
                if (inst->model_type == NULL_MODEL) {

                    if (strcmp(section_param, "TSP") != 0) {
                        printf("Only handles TSP. Aborted\n");
                        exit(EXIT_FAILURE);
                    }

                    /*TODO: add model types*/
                    inst->model_type = SYMMETRIC_TSP;
                }
                break;

            case DIMENSION:
                inst->num_nodes = atoi(section_param);
                inst->nodes = (node*) calloc(inst->num_nodes, sizeof(struct node_t));
                break;

            case EDGE_WEIGHT_TYPE:
                if(strcmp(section_param, "ATT") != 0) {
                    printf("Only handles ATT. Aborted\n");
                    /*exit(EXIT_FAILURE);*/ //  TODO: who cares
                }
                break;

            case NODE_COORD_SECTION:

                if(inst->num_nodes < 0) {
                    print_error("Define number of nodes before list of nodes. Aborted\n");
                }

                for (size_t i=1; i<=inst->num_nodes; i++) {
                    if ((getline(&line, &len, fp) == -1)) {
                        print_error("Reached end of file while reading nodes. Aborted\n");
                    }

                    size_t j,k;
                    size_t node_idx;
                    double x, y;

                    j = 0;
                    while (line[j] != ' ') j++;
                    line[j] = 0;
                    node_idx = atoi(line);

                    k = j+1;
                    while (line[j] != ' ') j++;
                    line[j] = 0;
                    x = atof(line+k);
                    y = atof(line+j+1);

                    if(node_idx != i) {
                        printf("Error in node indexing. Aborted\n");
                        exit(EXIT_FAILURE);
                    }

                    inst->nodes[i-1].x = x;
                    inst->nodes[i-1].y = y;
                }
                break;

            case END_OF_FILE:
                break;

            case TOUR_SECTION:
            case UNMANAGED_SECTION:
                printf("Warning: section %s unmanaged\n", section_name);
                break;
        }
        free(section_name);
        free(section_param);
	}

    free(line);
	fclose(fp);
    free(filename);
}

// TODO: merge with prec funct
solution parse_optimal_tour(instance inst) {

    char* filename;
    filename = (char*) calloc(100, sizeof(char));
    sprintf(filename, "../data/%s/%s.opt.tour", inst->model_name, inst->model_name);

	FILE *fp;
    fp = fopen(filename, "r");
	if (fp == NULL) {
        char* errormessage;
        errormessage = (char*) calloc(100, sizeof(char));
        sprintf(errormessage, "File %s not found. Aborted\n", filename);
        print_error(errormessage);
    }

    /* prepare tour node list */
    size_t* node_list;
    node_list = (size_t*) calloc(inst->num_nodes, sizeof(size_t));


	char* line;
    size_t len = 0;

	while((getline(&line, &len, fp) != -1)) {
		if (strlen(line) <= 1) continue;
        line[strlen(line)-1] = 0;  /* get rid of linefeed */

        /* retrive section and inject parameter */
        switch (hash2section(line)) {
            case TOUR_SECTION:
                /* read the nodes in the tour */
                for (size_t k=0; k<inst->num_nodes; k++) {
                    if ((getline(&line, &len, fp) == -1)) {
                        print_error("Reached end of file while reading nodes. Aborted\n");
                    }

                    node_list[k] = atoi(line)-1;
                }
                getline(&line, &len, fp); /* get rid of -1 */
                break;

            case NAME:
            case COMMENT:
            case TYPE:
            case DIMENSION:
            case EDGE_WEIGHT_TYPE:
            case NODE_COORD_SECTION:
            case UNMANAGED_SECTION:
                break;

            case END_OF_FILE:
                break;
        }
	}

    free(line);
	fclose(fp);
    free(filename);

	/* populate solution */
	solution sol = (solution) malloc(sizeof(struct solution_t));

	sol->optimality = OPTIMAL_TOUR;

	sol->num_edges = inst->num_nodes;
	sol->edges = (edge*) calloc(sol->num_edges, sizeof(struct edge_t));
    for (size_t k=0; k<sol->num_edges; k++) {
        sol->edges[k].i = node_list[k];
        sol->edges[k].j = node_list[(k+1) % sol->num_edges];
    }

    sol->zstar = 0.0;
    for (size_t k=0; k<sol->num_edges; k++) {
        sol->zstar += dist(sol->edges[k].i, sol->edges[k].j, inst);
    }

    return sol;
}

enum tsp_sections hash2section(char* section_name) {

    if (strcmp(section_name, "NAME") == 0)               return NAME;
    if (strcmp(section_name, "COMMENT") == 0)            return COMMENT;
    if (strcmp(section_name, "TYPE") == 0)               return TYPE;
    if (strcmp(section_name, "DIMENSION") == 0)          return DIMENSION;
    if (strcmp(section_name, "EDGE_WEIGHT_TYPE") == 0)   return EDGE_WEIGHT_TYPE;
    if (strcmp(section_name, "NODE_COORD_SECTION") == 0) return NODE_COORD_SECTION;
    if (strcmp(section_name, "TOUR_SECTION") == 0)       return TOUR_SECTION;
    if (strcmp(section_name, "EOF") == 0)                return END_OF_FILE;

    return UNMANAGED_SECTION;
}
