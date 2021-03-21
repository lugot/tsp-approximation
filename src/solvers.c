#include "tsp.h"
#include "utils.h"
#include "globals.h"
#include <cplex.h>
#include <stdlib.h>

// TODO:  retrive zstar from CPLEX
solution TSPopt(instance inst) {
	if (inst->model_type != SYMMETRIC_TSP) print_error("Need TSP instance");

	/* open CPLEX model */
	int error;
	/* enviorment: allocate memory */
	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP");

	/* populate enviorment with model data */
	build_tsp_model(inst, env, lp);

	if (CPXmipopt(env,lp)) print_error("CPXmipopt() error");

	/* store the optimal solution found by CPLEX */
	int ncols = CPXgetnumcols(env, lp);
	double* xstar = (double*) calloc(ncols, sizeof(double));
	if (CPXgetx(env, lp, xstar, 0, ncols-1)) print_error("CPXgetx() error");


	/* populate solution */
	solution sol = (solution) malloc(sizeof(struct solution_t));

	sol->optimality = OPTIMAL;
	/*sol->zstar  */
	sol->num_edges = inst->num_nodes;
	sol->edges = (edge*) calloc(sol->num_edges, sizeof(struct edge_t));

	/* index over selected edges */
	size_t k = 0;
	/* check for all edges if is selected */
	for (size_t i=0; i<inst->num_nodes; i++) for (size_t j=i+1; j<inst->num_nodes; j++) {
		if (xstar[xpos(i, j, inst)] > 0.5) { // TODO: add EPSILON check
			sol->edges[k].i = i;
			sol->edges[k].j = j;
			k++;
		}
	}

	if (k != sol->num_edges) print_error("Invalid number of edges selected. Aborted");

	free(xstar);
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);

	return sol;
}

solution optimal_tour(instance inst) {
	if (inst->model_type != OPTIMAL_TOUR) print_error("Need optimal tour instance");
	instance_tour tour = (instance_tour) inst;

	/* populate solution */
	solution sol = (solution) malloc(sizeof(struct solution_t));

	sol->optimality = OPTIMAL_TOUR;
	sol->num_edges = tour->num_nodes;
	sol->edges = (edge*) calloc(tour->num_nodes, sizeof(struct edge_t));

	size_t act;
	act = 0;

	size_t k = 0;
	while (tour->parent[act] != 0) {
		sol->edges[k].j = act;
		sol->edges[k].i = tour->parent[act];

		k++;
	}


	/*sol->zstar = 0.0;*/
	/*for (size_t k=0; k<sol->num_edges; k++) {*/
		/*sol->zstar += dist(sol->edges[k].i, sol->edges[k].j, inst);*/
	/*}*/

	return sol;
}
