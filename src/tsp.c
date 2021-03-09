#include "tsp.h"
#include <stdlib.h>
#include <cplex.h>

tsp_instance create_tsp_instance() {

    tsp_instance instance = (tsp_instance) malloc(sizeof(struct tsp_instance_t));

    instance->model_type = NULL_MODEL;
    instance->old_benders = 0;
    instance->randomseed = 0;
    instance->num_threads = 0;
    instance->timelimit = CPX_INFBOUND;
    instance->input_file = NULL;
    instance->node_file = NULL;
    instance->available_memory = 4096;
    instance->max_nodes = -1;
    instance->cutoff = CPX_INFBOUND;
    instance->costs_type = INTEGER;

    return instance;
}

/*TODO: comment function*/
/*TODO: expand print*/
void print_instance(tsp_instance instance) {
    printf("--- parameters ---");

    printf("model type: ");
    if (instance->model_type == NULL_MODEL)    printf("NULL_MODEL\n");
    if (instance->model_type == SYMMETRIC_TSP) printf("SYMMETRIC_TSP\n");
    if (instance->model_type == TOUR)          printf("TOUR\n");

    printf("input file: ");
    if (instance->input_file == NULL) printf("<NO INPUT FILE>\n");
    else                              printf("%s\n", instance->input_file);

    printf("node file: ");
    if (instance->node_file == NULL) printf("<NO NODE FILE>\n");
    else                             printf("%s\n", instance->node_file);

    printf("old benders: %d\n", instance->old_benders);
    printf("random seed: %d\n", instance->randomseed);
    printf("number of threads: %d\n", instance->num_threads);
    printf("time limit: %lf\n", instance->timelimit);
    printf("available memory: %d MB\n", instance->available_memory);

    printf("max branching nodes: ");
    if (instance->max_nodes == -1) printf("unlimited\n");
    else                           printf("%d\n", instance->max_nodes);

    printf("cutoff: %lf\n", instance->cutoff);

    printf("costs type: ");
    if (instance->costs_type == INTEGER) printf("INTEGER\n");
    if (instance->costs_type == DOUBLE)  printf("DOUBLE\n");

    printf("--- ---");
}
