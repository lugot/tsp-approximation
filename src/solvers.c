#include "tsp.h"
#include "utils.h"
#include "union_find.h"
#include <cplex.h>
#include <stdlib.h>
#include <sys/time.h>

int retrive_symmetric_solution(double* xstar, instance inst, solution sol);
int retrive_asymmetric_solution(double* xstar, instance inst, solution sol);

solution TSPopt(instance inst, enum model_types model_type) {
	if (inst->instance_type != TSP) print_error("need TSP instance");

	/* populate solution */
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
	CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit);
	CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0);
	CPXsetdblparam(env, CPX_PARAM_EPRHS, 1e-9);

	/* populate enviorment with model data */
	sol->build_time = build_tsp_model(inst, env, lp, model_type);


	struct timeval start, end;
	gettimeofday(&start, NULL);

	union_find uf = uf_create(inst->num_nodes);

	do {
		/* solve! */
		if (CPXmipopt(env,lp)) print_error("CPXmipopt() error");

		/* store the optimal solution found by CPLEX */
		int ncols = CPXgetnumcols(env, lp);
		double* xstar = (double*) calloc(ncols, sizeof(double));
		if (CPXgetx(env, lp, xstar, 0, ncols-1)) print_error("CPXgetx() error\n");

		int k;
		switch (model_type) {
			case SYMMETRIC:
<<<<<<< HEAD
				k = retrive_symmetric_solution(xstar, inst, sol);
				break;

			case SYMMETRIC_BENDERS:
				k = retrive_symmetric_solution(xstar, inst, sol);
				add_cool_subtour_elimination(inst, env, lp);
				k = retrive_symmetric_solution(xstar, inst, sol);
				if (VERBOSE) printf("[VERBOSE] Benders num sets: %d\n", uf->num_sets);
=======
			case SYMMETRIC_BENDERS:
				k = retrive_symmetric_solution(xstar, inst, sol, uf);
>>>>>>> parent of f2d3960 (added benders)
				break;
			case ASYMMETRIC_MTZ:
			case ASYMMETRIC_GG:
				k = retrive_asymmetric_solution(xstar, inst, sol);
				break;
			default:
				k = 0;
				break;
		}

		if (k != sol->num_edges) print_error("invalid number of edges\n");
		free(xstar);

	} while (model_type == SYMMETRIC_BENDERS && uf->num_sets != 1);


	gettimeofday(&end, NULL);
	long seconds = (end.tv_sec - start.tv_sec);
	long micros = ((seconds * 1000000) + end.tv_usec) - (start.tv_usec);
	sol->solve_time = micros / 1000.0;

	add_solution(inst, sol);

	// TODO: get this from CPLEX
	sol->zstar = zstar(inst, sol);

	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);

	return sol;
}

int retrive_symmetric_solution(double* xstar, instance inst, solution sol) {
	/* index over selected edges */
	int k = 0;
	/* check for all edges if is selected */
	for (int i=0; i<inst->num_nodes; i++) for (int j=i+1; j<inst->num_nodes; j++) {
		if (xstar[xpos(i, j, inst)] > 0.5) {
			sol->parent[i] = j;
			sol->edges[k++] = (edge) {i, j};
		}
	}

	return k;
}

int retrive_asymmetric_solution(double* xstar, instance inst, solution sol) {
	/* index over selected edges */
	int k = 0;
	/* check for all edges if is selected */
	for (int i=0; i<inst->num_nodes; i++) for (int j=0; j<inst->num_nodes; j++) {
		if (xstar[xxpos(i, j, inst)] > 0.5) {
			sol->parent[i] = j;
			sol->edges[k++] = (edge) {i, j};
		}
	}

	return k;
}
