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


void add_MTZ_subtour(instance inst, CPXENVptr env, CPXLPptr lp, solution sol)
{
	int u[inst->num_nodes - 1];

	for(int i=0; i<inst->num_nodes - 1; i++) 
	{
		u[i]=0;
	}

	int count = 1; 
	
	//node 0 is considered as starting node
	for(int k=0; k<sol->num_edges; k++) 
	{
			
		//Check if the subtour has been already filled
		if (sol->edges[k].i !=0 && u[sol->edges[k].i - 1] == 0)
		{


			u[sol->edges[k].i -1] = count;
			count++;			

			int start = sol->edges[k].i;
			int node = sol->edges[k].j;

			printf("Tour starts at %d\n", start);

			//Check the subtour
			while(node != start)
			{
				printf("... %d\n",node);

				if(sol->edges[node].i!=0)
				{
					u[node - 1] = count;
					count++;
				}
				node = sol->edges[node].j;
			}

		}
	}

	for(int i=0; i<inst->num_nodes - 1; i++) 
	{
		printf("%d\n",u[i]);
	}

	char continuous = 'C';
	/* lower and upper bound of the variable */
	double lb = 0.0;
	double ub = 1.0;


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
		if (CPXgetnumcols(env, lp)-1 != offset +i-1) {
			print_error(" wrong position for x var.s");
		}

	}


	/* right hand site scalar */
	double rhs = inst->num_nodes - 1.0;
	/* sense of the constraint:
	 * E: equality
	 * G: greater or equal
	 * L: less or equal */
	char sense = 'L';

	/* add one constraint at the time */
	for (size_t i=1; i<inst->num_nodes; i++) {
		for (size_t h = 1; h<inst->num_nodes; h++) {

			/* fetch last row position in current model (better: the new one)
			 * in CPLEX row == constraint */
			int lastrow = CPXgetnumrows(env, lp);

			/* write the constraint name inside CPLEX */
			sprintf(cname[0], "u(%ld,%ld)", h+1, i+1);

			/* add the new constraint (with coefficent zero, null lhs) */
			if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
				print_error("wrong CPXnewrows [degree]");
			}

			/* the graph is complete so we only skip self loops */
			if (i == h) continue;

			/* change the two u_s variables cofficient*/
			/* change the coefficent from 0 to 1.0 */
			if (CPXchgcoef(env, lp, lastrow, offset + i-1, u[i-1])) {
				print_error("wrong CPXchgcoef [degree]");
			}
			/* change the coefficent from 0 to 1.0 */
			if (CPXchgcoef(env, lp, lastrow, offset + h-1, u[h-1])) {
				print_error("wrong CPXchgcoef [degree]");
			}

			// uses big_M coefficient with M = num_nodes
			if (CPXchgcoef(env, lp, lastrow, xxpos(i, h, inst), inst -> num_nodes)) {
				print_error("wrong CPXchgcoef [degree]");
			}
			
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


