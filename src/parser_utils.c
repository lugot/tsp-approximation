#include "parser_utils.h"

#include "globals.h"
#include "tsp.h"

#include <unistd.h>
#include <getopt.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/*TODO: add this to a help*/
void print_usage() {
    printf("Usage: rectangle -file num -b num\n");
}

/*TODO: comment*/
/*TODO: add help*/
void parse_command_line(int argc, char** argv, tsp_instance instance) {

    static struct option long_options[] = {
        {"verbose",       no_argument,       NULL, 'v'},
        {"file",          required_argument, NULL, 'f'},
        {"time_limit",    required_argument, NULL, 'l'},
        {"model_type",    required_argument, NULL, 't'},
        {"old_benders",   required_argument, NULL, 'o'},
        {"seed",          required_argument, NULL, 's'},
        {"threads",       required_argument, NULL, 'T'},
        {"memory",        required_argument, NULL, 'M'},
        {"node_file",     required_argument, NULL, 'n'},
        {"max_nodes",     required_argument, NULL, 'N'},
        {"cutoff",        required_argument, NULL, 'c'},
        {"integer_costs", no_argument,       NULL, 'i'},
        {"help",          no_argument,       NULL, 'h'},
        {0,               0,                 NULL, 0  }
    };

    int long_index, opt;
    long_index = opt = 0;
    while ((opt = getopt_long(argc, argv,"f:l:t:o:s:T:M:n:N:C:c:",
                    long_options, &long_index)) != -1) {
        switch (opt) {
            case 'v' : VERBOSE = 1;
                break;
            case 'f' : strcpy(optarg, instance->input_file);
                break;
            case 'l' : instance->timelimit = atof(optarg);
                break;
            case 't' : instance->model_type = SYMMETRIC_TSP;
                /*if (strcmp(optarg, "TSP") == 0)  instance->model_type = SYMMETRIC_TSP;*/
                /*if (strcmp(optarg, "TOUR") == 0) instance->model_type = TOUR;*/
                break;
            case 'o' : instance->old_benders = atoi(optarg);
                break;
            case 's' : instance->randomseed = abs(atoi(optarg));
                break;
            case 'T' : instance->num_threads = atoi(optarg);
                break;
            case 'M' : instance->available_memory = atoi(optarg);
                break;
            case 'n' : strcpy(optarg, instance->node_file);
                break;
            case 'N' : instance->max_nodes = atoi(optarg);
                break;
            case 'C' : instance->cutoff = atof(optarg);
                break;
            case 'c' : instance->costs_type = INTEGER;
                break;
            case '?' :
                printf("Unknown option `-%c'.\n", optopt);
                print_usage();
                exit(EXIT_FAILURE);
        }
    }

    if (VERBOSE) print_instance(instance);
}



/*void read_input(instance *inst) // simplified CVRP parser, not all SECTIONs detected*/
/*{*/

    /*FILE *fin = fopen(inst->input_file, "r");*/
    /*if ( fin == NULL ) print_error(" input file not found!");*/

    /*inst->nnodes = -1;*/
    /*inst->depot = -1;*/
    /*inst->nveh = -1;*/

    /*char line[180];*/
    /*char *par_name;*/
    /*char *token1;*/
    /*char *token2;*/

    /*int active_section = 0; // =1 NODE_COORD_SECTION, =2 DEMAND_SECTION, =3 DEPOT_SECTION*/

    /*int do_print = ( VERBOSE >= 1000 );*/

    /*while ( fgets(line, sizeof(line), fin) != NULL )*/
    /*{*/
        /*if ( VERBOSE >= 2000 ) { printf("%s",line); fflush(NULL); }*/
        /*if ( strlen(line) <= 1 ) continue; // skip empty lines*/
        /*par_name = strtok(line, " :");*/
        /*if ( VERBOSE >= 3000 ) { printf("parameter \"%s\" ",par_name); fflush(NULL); }*/

        /*if ( strncmp(par_name, "NAME", 4) == 0 )*/
        /*{*/
            /*active_section = 0;*/
            /*continue;*/
        /*}*/

        /*if ( strncmp(par_name, "COMMENT", 7) == 0 )*/
        /*{*/
            /*active_section = 0;*/
            /*token1 = strtok(NULL, "");*/
            /*if ( VERBOSE >= 10 ) printf(" ... solving instance %s with model %d\n\n", token1, inst->model_type);*/
            /*continue;*/
        /*}*/

        /*if ( strncmp(par_name, "TYPE", 4) == 0 )*/
        /*{*/
            /*token1 = strtok(NULL, " :");*/
            /*if ( strncmp(token1, "CVRP",4) != 0 ) print_error(" format error:  only TYPE == CVRP implemented so far!!!!!!");*/
            /*active_section = 0;*/
            /*continue;*/
        /*}*/


        /*if ( strncmp(par_name, "DIMENSION", 9) == 0 )*/
        /*{*/
            /*if ( inst->nnodes >= 0 ) print_error(" repeated DIMENSION section in input file");*/
            /*token1 = strtok(NULL, " :");*/
            /*inst->nnodes = atoi(token1);*/
            /*if ( do_print ) printf(" ... nnodes %d\n", inst->nnodes);*/
            /*inst->demand = (double *) calloc(inst->nnodes, sizeof(double));*/
            /*inst->xcoord = (double *) calloc(inst->nnodes, sizeof(double));*/
            /*inst->ycoord = (double *) calloc(inst->nnodes, sizeof(double));*/
            /*active_section = 0;*/
            /*continue;*/
        /*}*/

        /*if ( strncmp(par_name, "CAPACITY", 8) == 0 )*/
        /*{*/
            /*token1 = strtok(NULL, " :");*/
            /*inst->capacity = atof(token1);*/
            /*if ( do_print ) printf(" ... vehicle capacity %lf\n", inst->capacity);*/
            /*active_section = 0;*/
            /*continue;*/
        /*}*/


        /*if ( strncmp(par_name, "VEHICLES", 8) == 0 )*/
        /*{*/
            /*token1 = strtok(NULL, " :");*/
            /*inst->nveh = atoi(token1);*/
            /*if ( do_print ) printf(" ... n. vehicles %d\n", inst->nveh);*/
            /*active_section = 0;*/
            /*continue;*/
        /*}*/


        /*if ( strncmp(par_name, "EDGE_WEIGHT_TYPE", 16) == 0 )*/
        /*{*/
            /*token1 = strtok(NULL, " :");*/
            /*if ( strncmp(token1, "EUC_2D", 6) != 0 ) print_error(" format error:  only EDGE_WEIGHT_TYPE == EUC_2D implemented so far!!!!!!");*/
            /*active_section = 0;*/
            /*continue;*/
        /*}*/

        /*if ( strncmp(par_name, "NODE_COORD_SECTION", 18) == 0 )*/
        /*{*/
            /*if ( inst->nnodes <= 0 ) print_error(" ... DIMENSION section should appear before NODE_COORD_SECTION section");*/
            /*active_section = 1;*/
            /*continue;*/
        /*}*/

        /*if ( strncmp(par_name, "DEMAND_SECTION", 14) == 0 )*/
        /*{*/
            /*if ( inst->nnodes <= 0 ) print_error(" ... DIMENSION section should appear before DEMAND_SECTION section");*/
            /*active_section = 2;*/
            /*continue;*/
        /*}*/

        /*if ( strncmp(par_name, "DEPOT_SECTION", 13) == 0 )*/
        /*{*/
            /*if ( inst->depot >= 0 ) print_error(" ... DEPOT_SECTION repeated??");*/
            /*active_section = 3;*/
            /*continue;*/
        /*}*/


        /*if ( strncmp(par_name, "EOF", 3) == 0 )*/
        /*{*/
            /*active_section = 0;*/
            /*break;*/
        /*}*/


        /*if ( active_section == 1 ) // within NODE_COORD_SECTION*/
        /*{*/
            /*int i = atoi(par_name) - 1;*/
            /*if ( i < 0 || i >= inst->nnodes ) print_error(" ... unknown node in NODE_COORD_SECTION section");*/
            /*token1 = strtok(NULL, " :,");*/
            /*token2 = strtok(NULL, " :,");*/
            /*inst->xcoord[i] = atof(token1);*/
            /*inst->ycoord[i] = atof(token2);*/
            /*if ( do_print ) printf(" ... node %4d at coordinates ( %15.7lf , %15.7lf )\n", i+1, inst->xcoord[i], inst->ycoord[i]);*/
            /*continue;*/
        /*}*/

        /*if ( active_section == 2 ) // within DEMAND_SECTION*/
        /*{*/
            /*int i = atoi(par_name) - 1;*/
            /*if ( i < 0 || i >= inst->nnodes ) print_error(" ... unknown node in NODE_COORD_SECTION section");*/
            /*token1 = strtok(NULL, " :,");*/
            /*inst->demand[i] = atof(token1);*/
            /*if ( do_print ) printf(" ... node %4d has demand %10.5lf\n", i+1, inst->demand[i]);*/
            /*continue;*/
        /*}*/

        /*if ( active_section == 3 ) // within DEPOT_SECTION*/
        /*{*/
            /*int i = atoi(par_name) - 1;*/
            /*if ( i < 0 || i >= inst->nnodes ) continue;*/
            /*if ( inst->depot >= 0 ) print_error(" ... multiple depots not supported in DEPOT_SECTION");*/
            /*inst->depot = i;*/
            /*if ( do_print ) printf(" ... depot node %d\n", inst->depot+1);*/
            /*continue;*/
        /*}*/

        /*printf(" final active section %d\n", active_section);*/
        /*print_error(" ... wrong format for the current simplified parser!!!!!!!!!");*/

    /*}*/

    /*fclose(fin);*/

/*}*/

