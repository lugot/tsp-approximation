#include "tsp.h"
#include "model_builder.h"
#include "globals.h"
#include "utils.h"
#include "union_find.h"
#include <cplex.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>

void retreive_symmetric_solution(double* xstar, solution sol);
void retreive_asymmetric_solution(double* xstar, solution sol);
/*int reachable(solution sol, int i, int j);*/
/*int reachable2(int* parent, int i, int j);*/
/*int fully_visit(solution sol);*/

int reachable(solution sol, int i, int j);
int visitable(solution sol);

/*static int CPXPUBLIC callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle);*/

solution TSPopt(instance inst, enum model_types model_type) {
	assert(inst->params != NULL && "no CPLEX params found");
	assert(inst->instance_type == TSP && "need TSP instance");

	/* create and populate solution */
	solution sol = (solution) calloc(1, sizeof(struct solution_t));
	sol->model_type = model_type;

	sol->num_edges = inst->num_nodes;
	sol->edges = (edge*) calloc(sol->num_edges, sizeof(struct edge_t));
	sol->parent = (int*) calloc(sol->num_edges, sizeof(int));

	sol->distance_time = compute_dist(inst);

	/* open CPLEX model */
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP");

	/* set params:
	 * - timelimit
	 * - EpInteger: CPLEX tollerance to declare a variable integer, important in bigM
	 * - EpRightHandSide: the less or equal satisfied up to this tollerance, usually 1e-5, with bigM 1e-9 */
	CPXsetdblparam(env, CPX_PARAM_TILIM, inst->params->timelimit);
	CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0);
	CPXsetdblparam(env, CPX_PARAM_EPRHS, 1e-9);

	/* populate enviorment with model data */
	sol->build_time = build_tsp_model(env, lp, inst, model_type);


	struct timeval start, end;
	gettimeofday(&start, NULL);

	/* solve! */
	if (CPXmipopt(env,lp)) print_error("CPXmipopt() error");

	/* store the optimal solution found by CPLEX */
	int ncols = CPXgetnumcols(env, lp);
	double* xstar = (double*) calloc(ncols, sizeof(double));
	if (CPXgetx(env, lp, xstar, 0, ncols-1)) print_error("CPXgetx() error\n");

	switch (model_type) {
		case SYMMETRIC:
			retreive_symmetric_solution(xstar, sol);
			break;

		case SYMMETRIC_BENDERS:
			retreive_symmetric_solution(xstar, sol);

			do {
				/* add the constraints */
				add_BENDERS_sec(env, lp, sol);

				/* solve! */
				if (CPXmipopt(env,lp)) print_error("CPXmipopt() error");

				/* store the optimal solution found by CPLEX */
				int ncols = CPXgetnumcols(env, lp);
				double* xstar = (double*) calloc(ncols, sizeof(double));
				if (CPXgetx(env, lp, xstar, 0, ncols-1)) print_error("CPXgetx() error\n");

				/* retrive the solution */
				retreive_symmetric_solution(xstar, sol);

			} while (!visitable(sol));

			/* save the complete model */
			char* filename;
			filename = (char*) calloc(100, sizeof(char));
			if (inst->model_folder == TSPLIB)    sprintf(filename, "../data_tsplib/%s/%s.benders.lp", inst->model_name, inst->model_name);
			if (inst->model_folder == GENERATED) sprintf(filename, "../data_generated/%s/%s.benders.lp", inst->model_name, inst->model_name);
			CPXwriteprob(env, lp, filename, NULL);
			free(filename);

			break;

		case ASYMMETRIC_MTZ:
		case ASYMMETRIC_GG:
			retreive_asymmetric_solution(xstar, sol);
			break;

		default:
			break;
	}

	free(xstar);

	gettimeofday(&end, NULL);
	long seconds = (end.tv_sec - start.tv_sec);
	long micros = ((seconds * 1000000) + end.tv_usec) - (start.tv_usec);
	sol->solve_time = micros / 1000.0;

	/* add the solution to the pool associated with it's instance */
	add_solution(inst, sol);

	sol->zstar = zstar(inst, sol); // TODO: add from cplex!

	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);

	return sol;
}

void retreive_symmetric_solution(double* xstar, solution sol) {
	/* set up parent as linked list */
	for (int i=0; i<sol->num_edges; i++) sol->parent[i] = i;

	/* index over selected edges */
	int k = 0;

	/* check for edges j>i if is selected */
	for (int i=0; i<sol->num_edges; i++) for (int j=i+1; j<sol->num_edges; j++) {

		if (xstar[xpos(i, j, sol->num_edges)] > 0.5) {
			sol->edges[k++] = (edge) {i, j};

			/* actually a shifted linked list */
			if (!reachable(sol, i, j)) swap(&sol->parent[i], &sol->parent[j]);
		}
	}

	assert(k == sol->num_edges && "not enought edges CPLEX solution");
}

void retreive_asymmetric_solution(double* xstar, solution sol) {
	/* set up parent as linked list */
	for (int i=0; i<sol->num_edges; i++) sol->parent[i] = i;

	/* index over selected edges */
	int k = 0;

	/* check for all edges if is selected */
	for (int i=0; i<sol->num_edges; i++) for (int j=0; j<sol->num_edges; j++) {

		if (xstar[xxpos(i, j, sol->num_edges)] > 0.5) {
			sol->edges[k++] = (edge) {i, j};

			/* actually a shifted linked list */
			if (!reachable(sol, i, j)) swap(&sol->parent[i], &sol->parent[j]);
		}
	}

	assert(k == sol->num_edges && "not enought edges CPLEX solution");
}

/* check if node j is reachable in the list by node i */
int reachable(solution sol, int i, int j) {
	int next = i;
	do {
		next = sol->parent[next];
		if (next == j) return 1;
	} while (sol->parent[next] != i);

	return 0;
}

/* check if all the nodes are visitated in a single tour*/
int visitable(solution sol) {
	int visits = 0;

	int next = 0;
	do {
		next = sol->parent[next];
		visits++;
	} while (sol->parent[next] != 0);

	return visits == sol->num_edges-1;
}
