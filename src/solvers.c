#include "tsp.h"
#include "utils.h"
#include <cplex.h>
#include <stdlib.h>
#include <sys/time.h>

// TODO:  retrive zstar from CPLEX
solution TSPopt(instance inst, enum model_types model_type) {
	if (inst->instance_type != TSP) print_error("need TSP instance");

	/* populate solution */
	solution sol = (solution) calloc(1, sizeof(struct solution_t));
	sol->model_type = model_type;
	sol->zstar = 0.0; // TODO: change

	sol->num_edges = inst->num_nodes;
	sol->edges = (edge*) calloc(sol->num_edges, sizeof(struct edge_t));
	sol->parent = (int*) calloc(sol->num_edges, sizeof(int));

	sol->distance_time = compute_dist(inst);


	/* open CPLEX model */
	int error;
	/* enviorment: allocate memory */
	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP");

	/* populate enviorment with model data */
	sol->build_time = build_tsp_model(inst, env, lp, model_type);

	struct timeval start, end;
    gettimeofday(&start, NULL);

	if (CPXmipopt(env,lp)) print_error("CPXmipopt() error");

    gettimeofday(&end, NULL);
    long seconds = (end.tv_sec - start.tv_sec);
    long micros = ((seconds * 1000000) + end.tv_usec) - (start.tv_usec);
	sol->solve_time = micros / 1000.0;

	/* store the optimal solution found by CPLEX */
	int ncols = CPXgetnumcols(env, lp);
	double* xstar = (double*) calloc(ncols, sizeof(double));
	if (CPXgetx(env, lp, xstar, 0, ncols-1)) print_error("CPXgetx() error");




	int (*cplex_pos)(int, int, instance) = model_type == SYMMETRIC ?
		xpos : xxpos;

	/* index over selected edges */
	int k = 0;
	/* check for all edges if is selected */
	for (int i=0; i<inst->num_nodes; i++) for (int j=0; j<inst->num_nodes; j++) {
		if (model_type == SYMMETRIC && j < i) continue;

		if (xstar[cplex_pos(i, j, inst)] > 0.5) {
			sol->parent[i] = j;
			sol->edges[k++] = (edge) {i, j};
		}
	}

	if (k != sol->num_edges) print_error("invalid number of edges");

	add_solution(inst, sol);

	free(xstar);
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);

	return sol;
}
