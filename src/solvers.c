#include "tsp.h"
#include "utils.h"
#include "globals.h"

#include <cplex.h>
#include <stdlib.h>

solution TSPopt(instance inst) {

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
	sol->num_edges = inst->num_nodes-1;
	sol->edges = (edge*) calloc(sol->num_edges, sizeof(struct edge_t));

	/* index over selected edges */
	size_t k = 0;
	/* check for all edges if is selected */
	for (size_t i=0; i<inst->num_nodes; i++) for (size_t j=i+1; j<inst->num_nodes; j++) {
		if (xstar[xpos(i, j, inst)] > 0.5) {
			sol->edges[k].i = i;
			sol->edges[k].j = j;
			k++;
		}
	}

	if (k-1 != sol->num_edges) print_error("Invalid number of edges selected. Aborted");

	free(xstar);
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);

	return sol;
}
