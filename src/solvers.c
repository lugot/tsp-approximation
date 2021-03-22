#include "tsp.h"
#include "utils.h"
#include <cplex.h>
#include <stdlib.h>

// TODO:  retrive zstar from CPLEX
solution TSPopt(instance inst, enum optimalities opt) {
	if (inst->type != TSP) print_error("need TSP instance");

	/* open CPLEX model */
	int error;
	/* enviorment: allocate memory */
	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP");

	/* populate enviorment with model data */
	build_tsp_model(inst, env, lp, opt);

	if (CPXmipopt(env,lp)) print_error("CPXmipopt() error");

	/* store the optimal solution found by CPLEX */
	int ncols = CPXgetnumcols(env, lp);
	double* xstar = (double*) calloc(ncols, sizeof(double));
	if (CPXgetx(env, lp, xstar, 0, ncols-1)) print_error("CPXgetx() error");


	/* populate solution */
	solution sol = (solution) calloc(1, sizeof(struct solution_t));
	sol->optimality = opt;
	sol->tsp_type = (opt == SYMMETRIC_SUBTOUR) ? SYMMETRIC : ASYMMETRIC;
	sol->zstar = 0.0; // TODO: change

	sol->num_edges = inst->num_nodes;
	sol->edges = (edge*) calloc(sol->num_edges, sizeof(struct edge_t));
	sol->parent = (int*) calloc(sol->num_edges, sizeof(int));

	int (*pos_checker)(int, int, instance) = opt == SYMMETRIC_SUBTOUR ?
		xpos : xxpos;

	/* index over selected edges */
	int k = 0;
	/* check for all edges if is selected */
	for (int i=0; i<inst->num_nodes; i++) for (int j=0; j<inst->num_nodes; j++) {
		if (opt == SYMMETRIC_SUBTOUR && j < i) continue;

		if (xstar[pos_checker(i, j, inst)] > 0.5) {
			sol->parent[i] = j;
			sol->edges[k++] = (edge) {i, j};
		}
	}

	if (opt == ASYMMETRIC_MTZ) {

		add_MTZ_subtour(inst, env, lp, sol);

		CPXwriteprob (env, lp, "myprob.mps", NULL);

		printf("Qui!\n");

		ncols = CPXgetnumcols(env, lp);
		double* xstar2 = (double*) calloc(ncols, sizeof(double));
		if (CPXgetx(env, lp, xstar2, 0, ncols-1)) print_error("CPXgetx() error");
	}

	if (k != sol->num_edges) print_error("invalid number of edges");

	add_solution(inst, sol);

	free(xstar);
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);

	return sol;
}
