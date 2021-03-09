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

	if (VERBOSE) printf("Reading input file %s\n", instance->input_file);

	FILE *fin = fopen(instance->input_file, "r");
	if (fin == NULL) {
        printf("File %s not found. Aborted\n", instance->input_file);
        exit(EXIT_FAILURE);
    }


	char* line = (char*) calloc(200, sizeof(char));
    char *model_name, *model_comment, *model_type;
	while(fgets(line, sizeof(line), fin) != NULL) {
		if (strlen(line) <= 1) continue;

		char* section_name = strtok(line, " : ");
        switch (hash2section(section_name)) {
            case NAME:
                model_name = strtok(NULL, " : ");
                assert(model_name != NULL);
                printf("nnnn |%s|\n", model_name);
                instance->model_name = (char*) malloc(strlen(model_name)*sizeof(char));
                strcpy(instance->model_name, "");
                break;

            case COMMENT:
                c_syntaxhell();
                model_comment = strtok(NULL, " : ");
                instance->model_comment = (char*) malloc(strlen(model_comment)*sizeof(char));
                strcpy(instance->model_comment, model_comment);
                break;

            case TYPE:
                /* keep the one set by user */
                if (instance->model_type == NULL_MODEL) {
                    model_type = strtok(NULL, " : ");

                    if (strncmp(model_type, "TSP", 3) != 0) {
                        printf("Only handles TSP. Aborted\n");
                        exit(EXIT_FAILURE);
                    }

                    /*TODO: add model types*/
                    instance->model_type = SYMMETRIC_TSP;
                }
                break;

            case DIMENSION:
                instance->num_nodes = atoi(strtok(NULL, " :"));

                instance->xcoord = (double*) calloc(instance->num_nodes, sizeof(double));
                instance->ycoord = (double*) calloc(instance->num_nodes, sizeof(double));
                break;

            case EDGE_WEIGHT_TYPE:
                if(strncmp(strtok(NULL, " : "), "ATT", 3) != 0) {
                    printf("Only handles ATT. Aborted\n");
                    exit(EXIT_FAILURE);
                }
                break;

            case NODE_COORD_SECTION:

                if(instance->num_nodes < 0) {
                    printf("Define number of nodes before list of nodes. Aborted\n");
                    exit(EXIT_FAILURE);
                }

                for (int i=1; i<=instance->num_nodes; i++) {
                    if (fgets(line, sizeof(line), fin) == NULL) {
                        printf("Reached end of file while reading nodes. Aborted\n");
                        exit(EXIT_FAILURE);
                    }

                    int node_index = atoi(strtok(line, " : "));
                    if(node_index != i) {
                        printf("Error in node indexing. Aborted\n");
                        exit(EXIT_FAILURE);
                    }

                    instance->xcoord[i-1] = atoi(strtok(NULL, " : "));
                    instance->ycoord[i-1] = atoi(strtok(NULL, " : "));
                }
                break;

            case END_OF_FILE:
                break;

            case UNMANAGED_SECTION:
                printf("Warning: section %s unmanaged\n", section_name);
                break;
        }
	}

    free(line);
	fclose(fin);
}
