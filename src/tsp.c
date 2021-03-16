#include "tsp.h"
#include "utils.h"
#include "globals.h"
#include <stdio.h>
#include <stdlib.h>
#include <cplex.h>

instance create_tsp_instance() {
    instance inst = (instance) calloc(1, sizeof(struct instance_t));

    /* initializing only the non-zero default parameters */
    inst->model_type = SYMMETRIC_TSP;
    inst->num_threads = -1;
    inst->timelimit = CPX_INFBOUND;
    inst->available_memory = 4096;
    inst->costs_type = REAL;

    return inst;
}

void build_tsp_model(instance inst, CPXENVptr env, CPXLPptr lp) {
	/*
	 * our model is a complete graph G=(V,E), |V|=n, |E|=m
	 * no self loop: m = n chooses 2 = n*(n-1)/2
	 *
	 * min sum e in E c_e * x_e
	 *
	 * sum e in delta(h) x_e = 2, forall h in V
	 * 0 <= x_e <= 1 integer, forall e in E
	 *
	 * typically i consider i < j
	 */


	/*double zero = 0.0;*/
	/* variable type:
	 * B: binary
	 * C: continuous */
	char binary = 'B';
	/* lower and upper bound of the variable */
	double lb = 0.0;
	double ub = 1.0;

	/* ColumnNAME: array of array used to inject variables in CPLEX */
	char** cname = (char**) calloc(1, sizeof(char *));
	cname[0] = (char*) calloc(100, sizeof(char));

	/* add a single binary variables x(i,j) for i < j at the time*/
	for (size_t i=0; i<inst->num_nodes; i++) for (size_t j=i+1; j<inst->num_nodes; j++) {

		/* write name of 1-indexed variable inside CPLEX */
		sprintf(cname[0], "x(%ld,%ld)", i+1, j+1);

		/* cost of the variable x(i, j): distance between nodes i and j */
		double obj = dist(i, j,inst);

		/* insert a single variable at the time so we can avoid passing arrays
		 * for the parameters
		 * in CPLEX variable == column */
		if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname)) {
			print_error("wrong CPXnewcols on x var.s");
		}

		/* test of xpos function */
		if (CPXgetnumcols(env, lp)-1 != xpos(i, j, inst)) {
			print_error(" wrong position for x var.s");
		}
	}


	/* right hand site scalar */
	double rhs = 2.0;
	/* sense of the constraint:
	 * E: equality
	 * G: greater or equal
	 * L: less or equal */
	char sense = 'E';

	/* add one constraint at the time */
	for (size_t h=0; h<inst->num_nodes; h++) {

		/* fetch last row position in current model (better: the new one)
		 * in CPLEX row == constraint */
		int lastrow = CPXgetnumrows(env, lp);

		/* write the constraint name inside CPLEX */
		sprintf(cname[0], "degree(%ld)", h+1);

		/* add the new constraint (with coefficent zero, null lhs) */
		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
			print_error("wrong CPXnewrows [degree]");
		}
		/* and change the coefficent if node i is adiacent to h */
		for (int i = 0; i<inst->num_nodes; i++) {
			/* the graph is complete so we only skip self loops */
			if (i == h) continue;

			/* change the coefficent from 0 to 1.0 */
			if (CPXchgcoef(env, lp, lastrow, xpos(i, h, inst), 1.0)) {
				print_error("wrong CPXchgcoef [degree]");
			}
		}
	}
	/* note: there are no subtour elimination constraints, relaxed model! */

	if (VERBOSE) CPXwriteprob(env, lp, "../data/relaxed_tsp_model.lp", NULL);

	free(cname[0]);
	free(cname);
}

void print_solution(solution sol) {
	printf("--- solution ---\n");
	printf("optimality: ");
	if (sol->optimality == OPTIMAL) printf("OPTIMAL\n");
	else 						    printf("APPROXIMATED\n");
	printf("num edges: %ld\n", sol->num_edges);
	for (size_t i=0; i<sol->num_edges; i++) {
		printf("x(%ld, %ld) = 1\n", sol->edges[i].i, sol->edges[i].j);
	}

    printf("--- ---\n\n");
}

void print_instance(instance inst) {
    printf("--- parameters ---\n");

    printf("model type: ");
    if (inst->model_type == NULL_MODEL)    printf("NULL_MODEL\n");
    if (inst->model_type == SYMMETRIC_TSP) printf("SYMMETRIC_TSP\n");
    printf("input file: ");
    if (inst->input_file == NULL) printf("<NO INPUT FILE>\n");
    else                              printf("%s\n", inst->input_file);
    printf("random seed: %d\n", inst->randomseed);
    printf("number of threads: %d\n", inst->num_threads);
    printf("time limit: %lf\n", inst->timelimit);
    printf("available memory: %d MB\n", inst->available_memory);
    printf("costs type: ");
    if (inst->costs_type == INTEGER) printf("INTEGER\n");
    if (inst->costs_type == REAL)    printf("REAL\n");
    if (inst->num_nodes == 0) {
        printf("skipping input data section, probably not parsed yet\n\n");
        return;
    }

    printf("--- input data ---\n");
    printf("model name: %s\n", inst->model_name);
    printf("model comment: %s\n", inst->model_comment);
    printf("model type: ");
    if (inst->model_type == NULL_MODEL)    printf("NULL_MODEL\n");
    if (inst->model_type == SYMMETRIC_TSP) printf("SYMMETRIC_TSP\n");
    printf("number of nodes: %d\n", inst->num_nodes);
    if (VERBOSE) {
        for (int i=0; i<inst->num_nodes; i++) {
            printf("%lf %lf\n", inst->nodes[i].x, inst->nodes[i].y);
        }
    }

    printf("--- ---\n\n");
}

void plot_solution_graphviz(instance inst, solution sol) {

	double box_size = 20.0;
	double max_coord = 0.0;
	for (size_t i=0; i<inst->num_nodes; i++) {
		max_coord = max(max_coord, inst->nodes[i].x);
		max_coord = max(max_coord, inst->nodes[i].y);
	}


	char* filename;
	filename = (char*) calloc(100, sizeof(char));
	sprintf(filename, "../data/%s.dot", inst->model_name);

	FILE* fp;
	fp = fopen(filename, "w");

	fprintf(fp, "graph %s {\n", inst->model_name);
	for (size_t i=0; i<inst->num_nodes; i++) {
		double plot_posx = inst->nodes[i].x / max_coord * box_size;
		double plot_posy = inst->nodes[i].y / max_coord * box_size;

		fprintf(fp, "\t%ld [ pos = \"%lf,%lf!\"]\n", i, plot_posx, plot_posy);
	}
	fprintf(fp, "\n");
	for (size_t k=0; k<sol->num_edges; k++) {
		fprintf(fp, "\t%ld -- %ld\n", sol->edges[k].i, sol->edges[k].j);
	}
	fprintf(fp, "}");


	fclose(fp);
	free(filename);
}
