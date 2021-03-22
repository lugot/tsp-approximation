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

	// TODO: add CPLEX memory management thingy


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

	if (VERBOSE) {
		char* filename;
		filename = (char*) calloc(100, sizeof(char));
		sprintf(filename, "../data/%s/%s.model.lp", inst->model_name, inst->model_name);

		CPXwriteprob(env, lp, filename, NULL);
		free(filename);
	}

	free(cname[0]);
	free(cname);
}



void build_asymmetric_tsp_model(instance inst, CPXENVptr env, CPXLPptr lp) {
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
	double self_ub = 0.0;

	/* ColumnNAME: array of array used to inject variables in CPLEX */
	char** cname = (char**) calloc(1, sizeof(char *));
	cname[0] = (char*) calloc(100, sizeof(char));

	/* add a single binary variables x(i,j) for i < j at the time*/
	for (size_t i=0; i<inst->num_nodes; i++) for (size_t j=0; j<inst->num_nodes; j++) {

		/* write name of 1-indexed variable inside CPLEX */
		sprintf(cname[0], "x(%ld,%ld)", i+1, j+1);

		/* cost of the variable x(i, j): distance between nodes i and j */
		double obj = dist(i, j,inst);

		/* insert a single variable at the time so we can avoid passing arrays
		 * for the parameters
		 * in CPLEX variable == column */
		if (i!=j)
		{
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname)) {
				print_error("wrong CPXnewcols on x var.s");
			}
		}
		else
			if (CPXnewcols(env, lp, 1, &obj, &lb, &self_ub, &binary, cname)) {
				print_error("wrong CPXnewcols on x var.s");
			}

		/* test of xxpos function */
		if (CPXgetnumcols(env, lp)-1 != xxpos(i, j, inst)) {
			print_error(" wrong position for x var.s");
		}
	}


	/* right hand site scalar */
	double rhs = 1.0;
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
			if (CPXchgcoef(env, lp, lastrow, xxpos(i, h, inst), 1.0)) {
				print_error("wrong CPXchgcoef [degree]");
			}
		}
	}
	for (size_t i=0; i<inst->num_nodes; i++) {

		/* fetch last row position in current model (better: the new one)
		 * in CPLEX row == constraint */
		int lastrow = CPXgetnumrows(env, lp);

		/* write the constraint name inside CPLEX */
		sprintf(cname[0], "degree(%ld)", inst->num_nodes + i+1);

		/* add the new constraint (with coefficent zero, null lhs) */
		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
			print_error("wrong CPXnewrows [degree]");
		}
		/* and change the coefficent if node i is adiacent to h */
		for (int h = 0; h<inst->num_nodes; h++) {
			/* the graph is complete so we only skip self loops */
			if (i == h) continue;

			/* change the coefficent from 0 to 1.0 */
			if (CPXchgcoef(env, lp, lastrow, xxpos(i, h, inst), 1.0)) {
				print_error("wrong CPXchgcoef [degree]");
			}
		}
	}
	
	add_MTZ_subtour(inst, env, lp);

	if (VERBOSE) {
		char* filename;
		filename = (char*) calloc(100, sizeof(char));
		sprintf(filename, "../data/%s/%s.model.lp", inst->model_name, inst->model_name);

		CPXwriteprob(env, lp, filename, NULL);
		free(filename);
	}

	free(cname[0]);
	free(cname);
}


void add_MTZ_subtour(instance inst, CPXENVptr env, CPXLPptr lp)
{

	char continuous = 'I';
	/* lower and upper bound of the variable */
	double lb = 1.0;
	double ub = inst->num_nodes-1;


	int offset = inst->num_nodes * inst->num_nodes;

	/* ColumnNAME: array of array used to inject variables in CPLEX */
	char** cname = (char**) calloc(1, sizeof(char *));
	cname[0] = (char*) calloc(100, sizeof(char));

	/* add a single binary variables x(i,j) for i < j at the time*/
	for (size_t i=1; i<inst->num_nodes; i++) {

		/* write name of 1-indexed variable inside CPLEX */
		sprintf(cname[0], "u(%ld)", i+1);

		/* cost of the variable x(i, j): distance between nodes i and j */
		double obj = 0;

		/* insert a single variable at the time so we can avoid passing arrays
		 * for the parameters
		 * in CPLEX variable == column */
		
		if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &continuous, cname)) {
			print_error("wrong CPXnewcols on x var.s");
		}
		
		/* test of xxpos function */
		if (CPXgetnumcols(env, lp)-1 != upos(i, inst)) {
			print_error(" wrong position for x var.s");
		}

	}

	int izero = 0;
	int index[3]; 
	double value[3];

	// add lazy constraints  1.0 * u_i - 1.0 * u_j + M * x_ij <= M - 1, for each arc (i,j) not touching node 0	
	double big_M = inst->num_nodes - 1.0;
	double rhs = big_M -1.0;
	char sense = 'L';
	int nnz = 3;
	for ( int i = 1; i < inst->num_nodes; i++ ) // excluding node 0
	{
		for ( int j = 1; j < inst->num_nodes; j++ ) // excluding node 0 
		{
			if ( i == j ) continue;
			sprintf(cname[0], "u-consistency for arc (%d,%d)", i+1, j+1);
			index[0] = upos(i,inst);	
			value[0] = 1.0;	
			index[1] = upos(j,inst);
			value[1] = -1.0;
			index[2] = xxpos(i,j,inst);
			value[2] = big_M;
			if ( CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index, value, cname) ) print_error("wrong CPXlazyconstraints() for u-consistency");
		}
	}
	
	// add lazy constraints 1.0 * x_ij + 1.0 * x_ji <= 1, for each arc (i,j) with i < j
	rhs = 1.0; 
	nnz = 2;
	for ( int i = 0; i < inst->num_nodes; i++ ) 
	{
		for ( int j = i+1; j < inst->num_nodes; j++ ) 
		{
			sprintf(cname[0], "SEC on node pair (%d,%d)", i+1, j+1);
			index[0] = xxpos(i,j,inst);
			value[0] = 1.0;
			index[1] = xxpos(j,i,inst);
			value[1] = 1.0;
			if ( CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index, value, cname) ) print_error("wrong CPXlazyconstraints on 2-node SECs");
		}
	}




}

void print_solution(solution sol) {
	printf("--- solution ---\n");
	printf("optimality: ");
	if (sol->optimality == OPTIMAL) printf("OPTIMAL\n");
	else 						    printf("APPROXIMATED\n");
	printf("num edges: %ld\n", sol->num_edges);
	for (size_t i=0; i<sol->num_edges; i++) {
		printf("x(%ld, %ld) = 1\n", (sol->edges[i].i+1), (sol->edges[i].j+1));
	}

    printf("--- ---\n\n");
}

void print_instance(instance inst) {
    printf("--- parameters ---\n");

    printf("model type: ");
    if (inst->model_type == NULL_MODEL)    printf("NULL_MODEL\n");
    if (inst->model_type == SYMMETRIC_TSP) printf("SYMMETRIC_TSP\n");
    if (inst->model_type == ASYMMETRIC_TSP) printf("ASYMMETRIC_TSP\n");
    printf("input file: ");
    if (inst->model_name == NULL) printf("<NO INPUT FILE>\n");
    else                          printf("%s\n", inst->model_name);
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
    if (inst->model_type == ASYMMETRIC_TSP) printf("ASYMMETRIC_TSP\n");
    if (inst->model_type == SYMMETRIC_TSP) printf("SYMMETRIC_TSP\n");
    printf("number of nodes: %ld\n", inst->num_nodes);
    if (VERBOSE) {
        for (int i=0; i<inst->num_nodes; i++) {
            printf("%lf %lf\n", inst->nodes[i].x, inst->nodes[i].y);
        }
    }

    printf("--- ---\n\n");
}

void plot_solutions_graphviz(instance inst, solution* sols, size_t num_sols) {
	double box_size = 20.0; // TODO: make proportional to the number of nodes
	double max_coord = 0.0;
	for (size_t i=0; i<inst->num_nodes; i++) {
		max_coord = max(max_coord, fabs(inst->nodes[i].x));
		max_coord = max(max_coord, fabs(inst->nodes[i].y));
	}

	char* filename;
	filename = (char*) calloc(100, sizeof(char));
	sprintf(filename, "../data/%s/%s.dot", inst->model_name, inst->model_name);

	FILE* fp;
	fp = fopen(filename, "w");

	fprintf(fp, "graph %s {\n", inst->model_name);
	fprintf(fp, "\tnode [shape=circle fillcolor=white]\n");
	for (size_t i=0; i<inst->num_nodes; i++) {
		double plot_posx = inst->nodes[i].x / max_coord * box_size;
		double plot_posy = inst->nodes[i].y / max_coord * box_size;

		fprintf(fp, "\t%ld [ pos = \"%lf,%lf!\"]\n", i, plot_posx, plot_posy);
	}
	fprintf(fp, "\n");

	for (size_t u=0; u<num_sols; u++) {
		for (size_t k=0; k<sols[u]->num_edges; k++) {
			fprintf(fp, "\t%ld -- %ld", sols[u]->edges[k].i, sols[u]->edges[k].j);

			if (sols[u]->optimality == OPTIMAL_TOUR) fprintf(fp, " [color = red]");
			fprintf(fp, "\n");
		}
	}
	fprintf(fp, "}");

	fclose(fp);
	free(filename);
}

void plot_solution_graphviz(instance inst, solution sol) {
	plot_solutions_graphviz(inst, (solution[]) {sol}, 1);
}


