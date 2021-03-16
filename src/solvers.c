#include "tsp.h"
#include "utils.h"
#include "globals.h"

#include <cplex.h>
#include <stdlib.h>

int TSPopt(instance inst) {

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

	inst->sol = (solution) malloc(sizeof(struct solution_t));
	inst->sol->optimality = OPTIMAL;
	inst->sol->num_edges = ncols;
	inst->sol->edges = (edge*) calloc(ncols, sizeof(struct edge_t));

	/*if (ncols != inst->num_nodes-1) print_error("Invalid number of solution. Aborted");*/
	// TODO: add zstar

	/* index over selected edges */
	size_t k = 0;
	/* check for all edges if is selected */
	for (size_t i=0; i<inst->num_nodes; i++) for (size_t j=i+1; j<inst->num_nodes; j++) {
		if (xstar[xpos(i, j, inst)] > 0.5) {
			inst->sol->edges[k].i = i;
			inst->sol->edges[k].j = j;
			k++;
		}
	}
	free(xstar);

	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);

	return EXIT_SUCCESS;

}
