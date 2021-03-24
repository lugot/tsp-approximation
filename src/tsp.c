#include "tsp.h"
#include "utils.h"
#include "globals.h"
#include "string.h"
#include <stdio.h>
#include <stdlib.h>
#include <cplex.h>
#include <sys/time.h>

void add_MTZ_subtour(instance inst, CPXENVptr env, CPXLPptr lp);
void add_GG_subtour(instance inst, CPXENVptr env, CPXLPptr lp);

instance create_tsp_instance() {
    instance inst = (instance) calloc(1, sizeof(struct instance_t));

    /* initializing only the non-zero default parameters */
    inst->instance_type = TSP;
    inst->num_threads = -1;
    inst->timelimit = CPX_INFBOUND;
    inst->available_memory = 4096;
    inst->costs_type = REAL;

    return inst;
}

instance duplicate_instance_parameters(instance inst) {
    instance newone = (instance) calloc(1, sizeof(struct instance_t));

    newone->model_name = (char*) calloc(strlen(inst->model_name), sizeof(char));
    strcpy(newone->model_name, inst->model_name);
    newone->model_comment = (char*) calloc(strlen(inst->model_comment), sizeof(char));
    strcpy(newone->model_comment, inst->model_comment);

    newone->instance_type = inst->instance_type;
    newone->randomseed = inst->randomseed;
    newone->num_threads = inst->num_threads;
    newone->timelimit = inst->timelimit;
    newone->available_memory = inst->available_memory;
    newone->costs_type = inst->costs_type;

    return newone;
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

double build_tsp_model(instance inst, CPXENVptr env, CPXLPptr lp, enum model_types model_type) {
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

	struct timeval start, end;
    gettimeofday(&start, NULL);

    /* double zero = 0.0; */
	/* variable type:
	 * B: binary
	 * C: continuous */
    char binary = 'B';
	/* lower and upper bound of the variable */
    double lb = 0.0;
    double ub = 1.0;
    int (*cplex_pos)(int, int, instance);

	/* ColumnNAME: array of array used to inject variables in CPLEX */
    char** cname = (char**) calloc(1, sizeof(char *));
    cname[0] = (char*) calloc(100, sizeof(char));

	/* add a single binary variables x(i,j) for i < j at the time */
    for (int i=0; i<inst->num_nodes; i++) for (int j=0; j<inst->num_nodes; j++) {

        if (model_type == SYMMETRIC && j <= i) continue;

		/* write name of 1-indexed variable inside CPLEX */
        sprintf(cname[0], "x(%d,%d)", i+1, j+1);

		/* cost of the variable x(i, j): distance between nodes i and j */
        double obj = dist(i, j,inst);

        switch (model_type) {
            case SYMMETRIC:
                ub = 1.0;
                cplex_pos = xpos;
                break;
            case ASYMMETRIC_MTZ:
            case ASYMMETRIC_GG:
                ub = (i == j) ? 0.0 : 1.0;
                cplex_pos = xxpos;
            case OPTIMAL_TOUR:
                break;
        }

		/* insert a single variable at the time so we can avoid passing arrays
		 * for the parameters
		 * in CPLEX variable == column */
        if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname)) {
            print_error("wrong CPXnewcols on x var.s");
        }

		/* test of cplex_position function */
        if (CPXgetnumcols(env, lp)-1 != cplex_pos(i, j, inst)) {
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

        switch (model_type) {
            case SYMMETRIC:
                rhs = 2.0;
                cplex_pos = xpos;
                break;
            case ASYMMETRIC_MTZ:
            case ASYMMETRIC_GG:
                rhs = 1.0;
                cplex_pos = xxpos;
            case OPTIMAL_TOUR:
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
            if (CPXchgcoef(env, lp, lastrow, cplex_pos(i, h, inst), 1.0)) {
                print_error("wrong CPXchgcoef [degree]");
            }
        }
    }
    if (model_type == ASYMMETRIC_MTZ || model_type == ASYMMETRIC_GG) {
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

        if (model_type == ASYMMETRIC_MTZ) add_MTZ_subtour(inst, env, lp);
        if (model_type == ASYMMETRIC_GG) add_GG_subtour(inst, env, lp);
    }

    if (VERBOSE) {
        // TODO: add modle type to filename
        char* filename;
        filename = (char*) calloc(100, sizeof(char));
        sprintf(filename, "../data/%s/%s.model.lp", inst->model_name, inst->model_name);

        CPXwriteprob(env, lp, filename, NULL);
        free(filename);
    }

    gettimeofday(&end, NULL);
    long seconds = (end.tv_sec - start.tv_sec);
    long micros = ((seconds * 1000000) + end.tv_usec) - (start.tv_usec);
	return micros / 1000.0;
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

void add_GG_subtour(instance inst, CPXENVptr env, CPXLPptr lp) {

	char integer = 'I';
	/* lower and upper bound of the variable */
	double lb = 0.0;
	double ub = inst->num_nodes-1;


	/* ColumnNAME: array of array used to inject variables in CPLEX */
	char** cname = (char**) calloc(1, sizeof(char *));
	cname[0] = (char*) calloc(100, sizeof(char));

	/* add a single binary variables x(i,j) for i < j at the time*/
	for (size_t i=0; i<inst->num_nodes; i++) for(size_t j=0; j<inst->num_nodes; j++){

		if(i==j) continue;
		/* write name of 1-indexed variable inside CPLEX */
		sprintf(cname[0], "y(%ld)(%ld)", i+1,j+1);

		/* cost of the variable x(i, j): distance between nodes i and j */
		double obj = 0;

		/* insert a single variable at the time so we can avoid passing arrays
		 * for the parameters
		 * in CPLEX variable == column */

		if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &integer, cname)) {
			print_error("wrong CPXnewcols on x var.s");
		}

		/* test of xxpos function */
		if (CPXgetnumcols(env, lp)-1 != ypos(i, j, inst)) {
			print_error(" wrong position for x var.s");
		}

	}

	double rhs = 0;
	char sense = 'L';

	for (size_t i=0; i<inst->num_nodes; i++)
		for (size_t h=0; h<inst->num_nodes; h++) {

			if (i == h) continue;

            /* fetch last row position in current model (better: the new one)
             * in CPLEX row == constraint */
            int lastrow = CPXgetnumrows(env, lp);

            /* write the constraint name inside CPLEX */
            sprintf(cname[0], "xy(%ld)(%ld)", i+1,h+1);

            /* add the new constraint (with coefficent zero, null lhs) */
            if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
                print_error("wrong CPXnewrows [degree]");
            }

            /* change the coefficent from 0 to 1.0 */
            if (CPXchgcoef(env, lp, lastrow, xxpos(i, h, inst), 1-inst->num_nodes)) {
                print_error("wrong CPXchgcoef [degree]");
            }

            /* change the coefficent from 0 to 1.0 */
            if (CPXchgcoef(env, lp, lastrow, ypos(i, h, inst), 1.0)) {
                print_error("wrong CPXchgcoef [degree]");
            }

    }


    rhs = inst->num_nodes - 1;
    sense = 'E';

    int lastrow = CPXgetnumrows(env, lp);

    /* write the constraint name inside CPLEX */
    sprintf(cname[0], "y1j");

    /* add the new constraint (with coefficent zero, null lhs) */
    if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
        print_error("wrong CPXnewrows [degree]");
        }

    for (size_t h=1; h<inst->num_nodes; h++) {

         	if (CPXchgcoef(env, lp, lastrow, ypos(0, h, inst), 1.0)) {
            	print_error("wrong CPXchgcoef [degree]");
        }
    }

    rhs = 1;
    sense = 'E';

    for (size_t j=1; j<inst->num_nodes; j++) {

            /* fetch last row position in current model (better: the new one)
             * in CPLEX row == constraint */
            int lastrow = CPXgetnumrows(env, lp);

            /* write the constraint name inside CPLEX */
            sprintf(cname[0], "y(%ld)", j+1);

            /* add the new constraint (with coefficent zero, null lhs) */
            if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
                print_error("wrong CPXnewrows [degree]");
            }

            /* and change the coefficent if node i is adiacent to h */
            for (int i = 0; i<inst->num_nodes; i++) {
                /* the graph is complete so we only skip self loops */
                if (i == j) continue;

                /* change the coefficent from 0 to 1.0 */
                if (CPXchgcoef(env, lp, lastrow, ypos(i, j, inst), 1.0)) {
                    print_error("wrong CPXchgcoef [degree]");
                }
            }

            /* and change the coefficent if node i is adiacent to h */
            for (int h = 0; h<inst->num_nodes; h++) {
                /* the graph is complete so we only skip self loops */
                if (h == j) continue;

                /* change the coefficent from 0 to 1.0 */
                if (CPXchgcoef(env, lp, lastrow, ypos(j, h, inst), -1.0)) {
                    print_error("wrong CPXchgcoef [degree]");
                }
            }
        }



}
