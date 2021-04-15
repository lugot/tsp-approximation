#include "solvers.h"
#include "tsp.h"
#include "model_builder.h"
#include "globals.h"
#include "utils.h"
#include "union_find.h"
#include <cplex.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>


solution TSPopt(instance inst, enum model_types model_type) {
	assert(inst != NULL);
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
	
	char* logfile;
	logfile = (char*) calloc(100, sizeof(char));
	if (inst->model_folder == TSPLIB)    sprintf(logfile, "../data_tsplib/%s/%slog.txt", inst->model_name, inst->model_name);
	if (inst->model_folder == GENERATED) sprintf(logfile, "../data_generated/%s/%slog.txt", inst->model_name, inst->model_name);
	
	CPXsetlogfilename(env, logfile, "w");
	
	CPXgetdettime(env,&sol->start);
	CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0);
	CPXsetdblparam(env, CPX_PARAM_EPRHS, 1e-9);

	/* populate enviorment with model data */
	sol->build_time = build_tsp_model(env, lp, inst, model_type);

	/* add callback if required */
	if(model_type == SYMMETRIC_BENDERS_CALLBACK) {
		CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE;
        if (CPXcallbacksetfunc(env, lp, contextid, add_BENDERS_sec_callback, sol)) {
            print_error("CPXcallbacksetfunc() error");
		}
	}


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
		case SYMMETRIC_BENDERS_CALLBACK:
			retreive_symmetric_solution(xstar, sol);
			break;

		case SYMMETRIC_BENDERS:
			retreive_symmetric_solution(xstar, sol);

			while (!visitable(sol)) {
				/* add the constraints */
				add_BENDERS_sec(env, lp, sol);

				/* some time has passed, updat the timelimit */
				gettimeofday(&end, NULL);
				CPXsetdblparam(env, CPX_PARAM_TILIM, inst->params->timelimit - (end.tv_sec - start.tv_sec));

				/* solve! */
				if (CPXmipopt(env,lp)) print_error("CPXmipopt() error");

				/* store the optimal solution found by CPLEX */
				int num_cols = CPXgetnumcols(env, lp);
				double* xstar = (double*) calloc(num_cols, sizeof(double));
				if (CPXgetx(env, lp, xstar, 0, num_cols-1)) print_error("CPXgetx() error\n");

				/* retrive the solution */
				retreive_symmetric_solution(xstar, sol);
			}


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

		case OPTIMAL_TOUR:
			assert(model_type != OPTIMAL_TOUR && "tried to solve an optimal tour instance");
			break;
	}

	free(xstar);

	CPXgetdettime(env,&sol->end);

	gettimeofday(&end, NULL);
	sol->solve_time = end.tv_sec - start.tv_sec; //sec
	//sol->solve_time = sol->end - sol->start; //ticks

	CPXgetobjval(env, lp, &sol->zstar);
	/*sol->zstar = zstar(inst, sol);*/

	/* add the solution to the pool associated with it's instance */
	add_solution(inst, sol);

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
			if (!reachable(sol, i, j) && !reachable(sol, j, i)) swap(&sol->parent[i], &sol->parent[j]);
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


void save_results(instance* insts, int num_instances) {
	assert(insts != NULL);
	assert(insts[0] != NULL);

	/* remove and create new fresh csv */
	remove("../results/results.csv");
	FILE* fp;
	fp = fopen("../results/results.csv", "w");
	assert(fp != NULL && "file not found while saving .csv");

	/* save the data */
	int num_models = insts[0]->num_solutions;
	fprintf(fp, "%d,", num_models);

	for (int i=0; i<num_models; i++) {
		enum model_types model_type = insts[0]->sols[i]->model_type;

		switch (model_type) {
			case SYMMETRIC:
			case OPTIMAL_TOUR:
				assert((model_type == SYMMETRIC || model_type == OPTIMAL_TOUR) &&
						"cannot retreive time from this kind of models");
				break;

			case ASYMMETRIC_MTZ:
				fprintf(fp, "MTZ");
				break;

			case ASYMMETRIC_GG:
				fprintf(fp, "GG");
				break;

			case SYMMETRIC_BENDERS:
				fprintf(fp, "BENDERS");
				break;

			case SYMMETRIC_BENDERS_CALLBACK: ;
				fprintf(fp, "BENDERS_CALLBACK");
				break;
		}

		if (i < num_models-1) fprintf(fp, ",");
		else                  fprintf(fp, "\n");
	}


	for (int i=0; i<num_instances; i++) {
		instance inst = insts[i];
		fprintf(fp, "%s,", inst->model_name);

		assert(inst->num_solutions == num_models && "missing some solutions");

		for (int j=0; j<num_models; j++) {
			assert(inst->sols[j]->model_type == insts[0]->sols[j]->model_type);

			if (j < num_models-1) fprintf(fp, "%lf,",  inst->sols[j]->solve_time);
			else                  fprintf(fp, "%lf\n", inst->sols[j]->solve_time);
		}

	}

	fclose(fp);

	/* generate the plot */
	//TODO: adjust timelimit
	system("python ../results/perprof.py -D , -T 3600 -S 2 -M 20 ../results/results.csv ../results/pp.pdf -P 'model comparisons'");
}