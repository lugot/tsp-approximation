#include "parser_utils.h"

#include "globals.h"
#include "tsp.h"

#include <unistd.h>
#include <getopt.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

enum tsp_sections {NAME, COMMENT, TYPE, DIMENSION, EDGE_WEIGHT_TYPE,
    NODE_COORD_SECTION, END_OF_FILE, UNMANAGED_SECTION};
enum tsp_sections hash2section(char* section_name);

/*TODO: add this to a help*/
void print_usage() {
    printf("Usage: rectangle -file num -b num\n");
}

/*TODO: comment*/
/*TODO: add help*/
void parse_command_line(int argc, char** argv, tsp_instance instance) {

    static struct option long_options[] = {
        {"verbose",       no_argument,       NULL, 'v'},
        {"input_file",    required_argument, NULL, 'f'},
        {"time_limit",    required_argument, NULL, 'l'},
        {"model_type",    required_argument, NULL, 't'},
        /*{"old_benders",   required_argument, NULL, 'o'},*/
        {"seed",          required_argument, NULL, 's'},
        {"threads",       required_argument, NULL, 'T'},
        {"memory",        required_argument, NULL, 'M'},
        /*{"node_file",     required_argument, NULL, 'n'},*/
        /*{"max_nodes",     required_argument, NULL, 'N'},*/
        /*{"cutoff",        required_argument, NULL, 'c'},*/
        {"integer_costs", no_argument,       NULL, 'i'},
        {"help",          no_argument,       NULL, 'h'},
        {0,               0,                 NULL, 0  }
    };

    int long_index, opt;
    long_index = opt = 0;
    while ((opt = getopt_long(argc, argv,":f:l:t:o:s:T:M:n:N:C:c:",
                    long_options, &long_index)) != -1) {
        switch (opt) {
            case 'v': VERBOSE = 1;
                break;

            case 'f':
                instance->input_file = (char*) malloc(strlen(optarg)*sizeof(char));
                strcpy(instance->input_file, optarg);
                break;

            case 'l': instance->timelimit = atof(optarg);
                break;

            case 't': instance->model_type = SYMMETRIC_TSP;
                /*if (strcmp(optarg, "TSP") == 0)  instance->model_type = SYMMETRIC_TSP;*/
                /*if (strcmp(optarg, "TOUR") == 0) instance->model_type = TOUR;*/
                break;

            /*case 'o': instance->old_benders = atoi(optarg);*/
                /*break;*/
            case 's': instance->randomseed = abs(atoi(optarg));
                break;

            case 'T': instance->num_threads = atoi(optarg);
                break;

            case 'M': instance->available_memory = atoi(optarg);
                break;

            /*case 'n':*/
                /*instance->node_file = (char*) malloc(strlen(optarg)*sizeof(char));*/
                /*strcpy(instance->node_file, optarg);*/
                /*break;*/
            /*case 'N': instance->max_nodes = atoi(optarg);*/
                /*break;*/
            /*case 'C': instance->cutoff = atof(optarg);*/
                /*break;*/
            case 'c': instance->costs_type = INTEGER;
                break;

            case '?':
                printf("Unknown option `-%c'.\n", optopt);
                print_usage();
                exit(EXIT_FAILURE);
        }
    }
}


enum tsp_sections hash2section(char* section_name) {
    if (strcmp(section_name, "NAME") == 0)               return NAME;
    if (strcmp(section_name, "COMMENT") == 0)            return COMMENT;
    if (strcmp(section_name, "TYPE") == 0)               return TYPE;
    if (strcmp(section_name, "DIMENSION") == 0)          return DIMENSION;
    if (strcmp(section_name, "EDGE_WEIGHT_TYPE") == 0)   return EDGE_WEIGHT_TYPE;
    if (strcmp(section_name, "NODE_COORD_SECTION") == 0) return NODE_COORD_SECTION;
    if (strcmp(section_name, "EOF") == 0)                return END_OF_FILE;

    return UNMANAGED_SECTION;
}

/*TODO: comment function*/
void parse_input(tsp_instance instance) {

	if (VERBOSE) printf("[VERBOSE] reading input file %s\n", instance->input_file);

	FILE *fin = fopen(instance->input_file, "r");
	if (fin == NULL) {
        printf("File %s not found. Aborted\n", instance->input_file);
        exit(EXIT_FAILURE);
    }


	char* line;
    char *section_name, *section_param;
    size_t len = 0;

	while((getline(&line, &len, fin) != -1)) {
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
                instance->model_name = (char*) malloc(strlen(section_param)*sizeof(char));
                strcpy(instance->model_name, section_param);
                break;

            case COMMENT:
                instance->model_comment = (char*) malloc(strlen(section_param)*sizeof(char));
                strcpy(instance->model_comment, section_param);
                break;

            case TYPE:
                /* keep the one set by user */
                if (instance->model_type == NULL_MODEL) {

                    if (strcmp(section_param, "TSP") != 0) {
                        printf("Only handles TSP. Aborted\n");
                        exit(EXIT_FAILURE);
                    }

                    /*TODO: add model types*/
                    instance->model_type = SYMMETRIC_TSP;
                }
                break;

            case DIMENSION:
                instance->num_nodes = atoi(section_param);

                instance->xcoord = (double*) calloc(instance->num_nodes, sizeof(double));
                instance->ycoord = (double*) calloc(instance->num_nodes, sizeof(double));
                break;

            case EDGE_WEIGHT_TYPE:
                if(strcmp(section_param, "ATT") != 0) {
                    printf("Only handles ATT. Aborted\n");
                    exit(EXIT_FAILURE);
                }
                break;

            case NODE_COORD_SECTION:

                if(instance->num_nodes < 0) {
                    printf("Define number of nodes before list of nodes. Aborted\n");
                    exit(EXIT_FAILURE);
                }

                for (size_t i=1; i<=instance->num_nodes; i++) {
                    if ((getline(&line, &len, fin) == -1)) {
                        printf("Reached end of file while reading nodes. Aborted\n");
                        exit(EXIT_FAILURE);
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

                    instance->xcoord[i-1] = x;
                    instance->ycoord[i-1] = y;
                }
                break;

            case END_OF_FILE:
                break;

            case UNMANAGED_SECTION:
                printf("Warning: section %s unmanaged\n", section_name);
                break;
        }
        free(section_name);
        free(section_param);
	}

    free(line);
	fclose(fin);
}
