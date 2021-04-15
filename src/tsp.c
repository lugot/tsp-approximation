#include "tsp.h"
#include "string.h"
#include "utils.h"
#include <time.h>
#include <assert.h>

instance create_empty_instance() {
    instance inst = (instance) calloc(1, sizeof(struct instance_t));
    inst->params = (cplex_params) calloc(1, sizeof(struct cplex_params_t));

    /* initializing only the non-zero default parameters */
    inst->params->num_threads = -1;
    inst->params->timelimit = CPX_INFBOUND;
    inst->params->available_memory = 4096;
    inst->params->cost = REAL;

    return inst;
}

instance create_instance(cplex_params params) {
    instance inst = (instance) calloc(1, sizeof(struct instance_t));

    /* initializing with passed parameters */
    inst->params = params;

    return inst;
}

void add_params(instance inst, cplex_params params) {
    /* should be non NULL by default constructor */
    memcpy(inst->params, params, sizeof(struct cplex_params_t));
}

instance generate_random_instance(int id, int num_nodes, int box_size) {
    instance inst = create_empty_instance();
    srand(id);

    char* buf;
    buf = (char*) calloc(100, sizeof(char));

    sprintf(buf, "random%d", id);
    inst->model_name = (char*) calloc(strlen(buf), sizeof(char));
    strcpy(inst->model_name, buf);

    sprintf(buf, "random generated with id %d", id);
    inst->model_comment = (char*) calloc(strlen(buf), sizeof(char));
    strcpy(inst->model_comment, buf);

    inst->weight_type = ATT;
    inst->num_nodes = num_nodes;

    inst->nodes = (node*) calloc(num_nodes, sizeof(struct node_t));
    for (int i=0; i<inst->num_nodes; i++) {

        inst->nodes[i].x = (double) rand() / ((double) RAND_MAX/box_size);
        inst->nodes[i].y = (double) rand() / ((double) RAND_MAX/box_size);

        /* TODO: investigate
         *uint64_t r53 = ((uint64_t)(rand()) << 21) ^ (rand() >> 2);
         *return (double)r53 / 9007199254740991.0; // 2^53 - 1
         */
    }

    return inst;
}

instance* generate_random_instances(int num_instances, int num_nodes, int box_size) {
    instance* insts = (instance*) calloc(num_instances, sizeof(struct instance_t));

    //TODO: do it better
    srand(time(NULL));

    for (int i=0; i<num_instances; i++) {
        int rand_id = rand(); /* between 0 and RAND_MAX (2^31 or smht here) */

        insts[i] = generate_random_instance(rand_id, num_nodes, box_size);
    }

    return insts;
}

// TODO: test
instance clone_instance(instance inst) {
    instance newone = (instance) calloc(1, sizeof(struct instance_t));

    newone->model_name = (char*) calloc(strlen(inst->model_name), sizeof(char));
    strcpy(newone->model_name, inst->model_name);
    newone->model_comment = (char*) calloc(strlen(inst->model_comment), sizeof(char));
    strcpy(newone->model_comment, inst->model_comment);

    newone->params = (cplex_params) calloc(1, sizeof(struct cplex_params_t));
    memcpy(newone->params, inst->params, sizeof(struct cplex_params_t));

    return newone;
}

void free_instance(instance inst) {
    free(inst->model_name);
    free(inst->model_comment);

    free(inst->nodes);
    for (int i=0; i<inst->num_nodes; i++) free(inst->adjmatrix[i]);
    free(inst->adjmatrix);

    // TODO: free sols
    /*for (int i=0; i<inst->num_nodes; i++) free(inst->adjmatrix[i]);*/
    free(inst->sols);
}

void save_instance(instance inst) {
	char* filename;
	filename = (char*) calloc(100, sizeof(char));
	sprintf(filename, "../data_generated/%s/%s.tsp", inst->model_name, inst->model_name);
}


void add_solution(instance inst, solution sol) {
    inst->sols = (solution*) realloc(inst->sols, inst->num_solutions+1);
    inst->sols[inst->num_solutions++] = sol;
    sol->inst = inst;
}

/* printers */
void print_instance(instance inst, int print_data) {
	printf("\ninfos:\nmodel: %s\n", inst->model_name);
	printf("comment: %s\n", inst->model_comment);
	printf("- model type: ");
	switch (inst->instance_type) {
		case TSP:
			printf("TSP\n");
			break;
		case TOUR:
			printf("(solution only)\n");
			break;
		default:
			printf("\n");
			break;
	}

    print_cplex_params(inst->params);

	printf("input data:\n");
	printf("- weight type: ");
	switch (inst->weight_type) {
		case ATT:
			printf("ATT\n");
			break;
		case EUC_2D:
			printf("EUC_2D\n");
			break;
		case GEO:
			printf("GEO\n");
			break;
		case EXPLICIT:
			printf("(matrix explicit)\n");
			break;
		default:
			break;
	}
	printf("- number of nodes: %d\n", inst->num_nodes);
	if (print_data) {
		printf("- weights:\n");
		if (inst->adjmatrix != NULL) {
			/* compute numof figures for spacing */
			int column_width = 0;
			for (int i=0; i<inst->num_nodes; i++) for (int j=0; j<=i; j++) {
				column_width = max(column_width, 1 + (int) log10(inst->adjmatrix[i][j]));
			}
			column_width = 2 + max(column_width, 1 + (int) log10(inst->num_nodes));

			/* figures.ab so max 1 decimal figures */
			char* buffer = (char*) calloc(column_width, sizeof(char));

			for (int i=0; i<inst->num_nodes; i++) {
				sprintf(buffer, "%d | ", i+1);
				printf("%*s", column_width+1, buffer);

				for (int j=0; j<=i; j++) {
					sprintf(buffer, "%.1lf ", inst->adjmatrix[i][j]);
					printf("%*s", column_width+1, buffer);
				}
				printf("\n");
			}
			free(buffer);
		}
		printf("- nodes:\n");
		if (inst->nodes != NULL) {
			int column_width = 0;
			for (int i=0; i<inst->num_nodes; i++) {
				column_width = max(column_width, 1 + (int) log10(inst->nodes[i].x));
				column_width = max(column_width, 1 + (int) log10(inst->nodes[i].y));
			}
			column_width = 2 + max(column_width, 1 + (int) log10(inst->num_nodes));

			/* figures.ab so max 1 decimal figures */
			char* buffer = (char*) calloc(column_width, sizeof(char));

			for (int i=0; i<inst->num_nodes; i++) {
				sprintf(buffer, "%d | ", i+1);
				printf("%*s", column_width+1, buffer);

				sprintf(buffer, "%.1lf ", inst->nodes[i].x);
				printf("%*s", column_width+1, buffer);

				sprintf(buffer, "%.1lf ", inst->nodes[i].y);
				printf("%*s\n", column_width+1, buffer);
			}
			free(buffer);
		}
	}

	printf("solutions:\n");
	printf("- num of solutions: %d\n", inst->num_solutions);
	for (int i=0; i<inst->num_solutions; i++) {
		printf("solution %d:\n", i+1);
		print_solution(inst->sols[i], print_data);
	}

	printf("--- ---\n\n");
}

void print_cplex_params(cplex_params params) {
	printf("cplex params:\n");
    if (params == NULL) return;
	printf("- random seed: %d\n", params->randomseed);
	printf("- number of threads: %d\n", params->num_threads);
	printf("- time limit: %lf\n", params->timelimit);
	printf("- available memory: %d MB\n", params->available_memory);
	printf("- costs type: ");
	switch (params->cost) {
		case REAL:
			printf("real\n");
			break;
		case INTEGER:
			printf("integer\n");
			break;
	}
}

void print_solution(solution sol, int print_data) {
	printf("- optimality: ");
	switch (sol->model_type) {
		case OPTIMAL_TOUR:
			printf("optimal tour from file\n");
			break;
		case SYMMETRIC:
			printf("symmetric, no subtour elimination\n");
			break;
		case ASYMMETRIC_MTZ:
			printf("asymmetric, mtz subtour elimination\n");
			break;
		case ASYMMETRIC_GG:
			printf("asymmetric, gg subtour elimination\n");
			break;
		case ASYMMETRIC_PROF_GG:
			printf("asymmetric, gg subtour elimination\n");
			break;
		default:
			printf("\n");
			break;
	}
	printf("- zstar: %lf\n", sol->zstar);
	printf("- num edges: %d\n", sol->num_edges);
	if (print_data) {
		printf("- edges:\n");
		if (sol->edges != NULL) {
			int column_width = 2 + (int) log10(sol->num_edges);
			char* buffer = (char*) calloc(column_width, sizeof(char));

			for (int i=0; i<sol->num_edges; i++) {
				sprintf(buffer, "%d | ", i+1);
				printf("%*s", column_width+2, buffer);

				sprintf(buffer, "%d ", sol->edges[i].i+1);
				printf("%*s", column_width, buffer);

				sprintf(buffer, "%d ", sol->edges[i].j+1);
				printf("%*s\n", column_width, buffer);
			}
			free(buffer);
		}
        printf("- parent:\n");
        if (sol->parent != NULL) {
            int column_width = 1 + (int) log10(sol->num_edges);
            char* buffer = (char*) calloc(column_width, sizeof(char));

            for (int i=0; i<sol->num_edges; i++) {
                sprintf(buffer, "%d", i+1);
                printf("%*s", column_width+1, buffer);
            }
            printf("\n");
            for (int i=0; i<sol->num_edges; i++) {
                sprintf(buffer, "%d", sol->parent[i]+1);
                printf("%*s", column_width+1, buffer);
            }
            printf("\n");
            free(buffer);
        }
	}
	printf("- distance time: %lf ms\n", sol->distance_time);
	printf("- build time: %lf ms\n", sol->build_time);
	printf("- solve time: %lf s\n", sol->solve_time / 1000.0);
	printf("- solve time(ticks): %lf ticks\n", sol->end-sol->start);
}

void plot_solution_graphviz(solution sol) {
	plot_solutions_graphviz((solution[]) {sol}, 1);
}

void plot_solutions_graphviz(solution* sols, int num_sols) {
	double box_size = 20.0; // TODO: make proportional to the number of nodes
	double max_coord = 0.0;

	instance inst = sols[0]->inst;
    assert(inst != NULL && "no instance associated with solution");

	for (int i=0; i<inst->num_nodes && inst->nodes != NULL; i++) {
		max_coord = max(max_coord, fabs(inst->nodes[i].x));
		max_coord = max(max_coord, fabs(inst->nodes[i].y));
	}

	char* filename;
	filename = (char*) calloc(100, sizeof(char));
	if (inst->model_folder == TSPLIB)    sprintf(filename, "../data_tsplib/%s/%s.dot", inst->model_name, inst->model_name);
	if (inst->model_folder == GENERATED) sprintf(filename, "../data_generated/%s/%s.dot", inst->model_name, inst->model_name);

	FILE* fp;
	fp = fopen(filename, "w");
	assert(fp != NULL && "file not found while saving .dot");

	fprintf(fp, "graph %s {\n", inst->model_name);
	fprintf(fp, "\tnode [shape=circle fillcolor=white]\n");
	for (int i=0; i<inst->num_nodes && inst->nodes != NULL; i++) {
		double plot_posx = inst->nodes[i].x / max_coord * box_size;
		double plot_posy = inst->nodes[i].y / max_coord * box_size;

		fprintf(fp, "\t%d [ pos = \"%lf,%lf!\"]\n", i+1, plot_posx, plot_posy);
	}
	fprintf(fp, "\n");

	for (int u=0; u<num_sols; u++) {
		for (int k=0; k<sols[u]->num_edges; k++) {
			fprintf(fp, "\t%d -- %d", sols[u]->edges[k].i+1, sols[u]->edges[k].j+1);

			if (sols[u]->model_type == OPTIMAL_TOUR) fprintf(fp, " [color = red]");
			fprintf(fp, "\n");
		}
	}
	fprintf(fp, "}");

	fclose(fp);
	free(filename);
}




