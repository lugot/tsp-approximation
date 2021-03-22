#include "tsp.h"
#include "utils.h"
#include "globals.h"
#include "string.h"
#include <stdio.h>
#include <stdlib.h>
#include <cplex.h>

instance create_tsp_instance() {
    instance inst = (instance) calloc(1, sizeof(struct instance_t));

    /* initializing only the non-zero default parameters */
    inst->type = TSP;
    inst->num_threads = -1;
    inst->timelimit = CPX_INFBOUND;
    inst->available_memory = 4096;
    inst->costs_type = REAL;

    return inst;
}

instance duplicate_instance(instance inst) {
    instance newone = (instance) calloc(1, sizeof(struct instance_t));

    newone->model_name = (char*) calloc(strlen(inst->model_name), sizeof(char));
    strcpy(newone->model_name, inst->model_name);
    newone->model_comment = (char*) calloc(strlen(inst->model_comment), sizeof(char));
    strcpy(newone->model_comment, inst->model_comment);

    // TODO: add duplication
    /*newone->type = inst->type;*/
    /*newone->randomseed = inst->type;*/
    /*newone->num_threads = inst->type;*/
    /*newone->timelimit = inst->type;*/
    /*newone->= inst->type;*/


    /*[> initializing only the non-zero default parameters <]*/
    /*inst->type = TSP;*/
    /*inst->num_threads = -1;*/
    /*inst->timelimit = CPX_INFBOUND;*/
    /*inst->available_memory = 4096;*/
    /*inst->costs_type = REAL;*/

    return inst;
}
void free_intance(instance inst) {
    free(inst->model_name);
    free(inst->model_comment);

    free(inst->nodes);
    for (int i=0; i<inst->num_nodes; i++) free(inst->adjmatrix[i]);
    free(inst->adjmatrix);

    // TODO: free sols
    /*for (int i=0; i<inst->num_nodes; i++) free(inst->adjmatrix[i]);*/
    free(inst->sols);
}

void add_solution(instance inst, solution sol) {
    inst->sols = (solution*) realloc(inst->sols, inst->num_solutions+1);
    inst->sols[inst->num_solutions++] = sol;
    sol->inst = inst;
}

void build_tsp_model(instance inst, CPXENVptr env, CPXLPptr lp, enum optimalities opt) {
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

    /* double zero = 0.0; */
	/* variable type:
	 * B: binary
	 * C: continuous */
    char binary = 'B';
	/* lower and upper bound of the variable */
    double lb = 0.0;
    double ub = 1.0;
    int (*pos_checker)(int, int, instance);

	/* ColumnNAME: array of array used to inject variables in CPLEX */
    char** cname = (char**) calloc(1, sizeof(char *));
    cname[0] = (char*) calloc(100, sizeof(char));

	/* add a single binary variables x(i,j) for i < j at the time */
    for (int i=0; i<inst->num_nodes; i++) for (int j=0; j<inst->num_nodes; j++) {

        if (opt == SYMMETRIC_SUBTOUR && j <= i) continue;

		/* write name of 1-indexed variable inside CPLEX */
        sprintf(cname[0], "x(%d,%d)", i+1, j+1);

		/* cost of the variable x(i, j): distance between nodes i and j */
        double obj = dist(i, j,inst);

        switch (opt) {
            case SYMMETRIC_SUBTOUR:
                ub = 1.0;
                pos_checker = xpos;
                break;
            case ASYMMETRIC_MTZ:
                ub = (i == j) ? 0.0 : 1.0;
                pos_checker = xxpos;
            default:
                break;
        }

		/* insert a single variable at the time so we can avoid passing arrays
		 * for the parameters
		 * in CPLEX variable == column */
        if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname)) {
            print_error("wrong CPXnewcols on x var.s");
        }

		/* test of xpos function */
        if (CPXgetnumcols(env, lp)-1 != pos_checker(i, j, inst)) {
            print_error("wrong position for x var.s");
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
    for (int h=0; h<inst->num_nodes; h++) {

		/* fetch last row position in current model (better: the new one)
		 * in CPLEX row == constraint */
        int lastrow = CPXgetnumrows(env, lp);

		/* write the constraint name inside CPLEX */
        sprintf(cname[0], "degree(%d)", h+1);

        switch (opt) {
            case SYMMETRIC_SUBTOUR:
                rhs = 2.0;
                pos_checker = xpos;
                break;
            case ASYMMETRIC_MTZ:
                rhs = 1.0;
                pos_checker = xxpos;
            default:
                break;
        }


		/* add the new constraint (with coefficent zero, null lhs) */
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
            print_error("wrong CPXnewrows [degree]");
        }
		/* and change the coefficent if node i is adiacent to h */
        for (int i = 0; i<inst->num_nodes; i++) {
            /* the graph is complete so we only skip self loops */
            if (i == h) continue;

            /* change the coefficent from 0 to 1.0 */
            if (CPXchgcoef(env, lp, lastrow, pos_checker(i, h, inst), 1.0)) {
                print_error("wrong CPXchgcoef [degree]");
            }
        }
    }
    if (opt == ASYMMETRIC_MTZ) {
        /* note: there are no subtour elimination constraints, relaxed model! */
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
    }

    if (VERBOSE) {
        // TODO: add modle type to filename
        char* filename;
        filename = (char*) calloc(100, sizeof(char));
        sprintf(filename, "../data/%s/%s.model.lp", inst->model_name, inst->model_name);

        CPXwriteprob(env, lp, filename, NULL);
        free(filename);
    }

    free(cname[0]);
    free(cname);
}

void add_MTZ_subtour(instance inst, CPXENVptr env, CPXLPptr lp) {

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

