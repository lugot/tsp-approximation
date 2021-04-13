#include "solvers.h"
#include "tsp.h"
#include "union_find.h"
#include "utils.h"
#include "globals.h"
#include "solvers.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <cplex.h>
#include <sys/time.h>
#include <time.h>
#include <assert.h>

void add_symmetric_variables(CPXENVptr env, CPXLPptr lp, instance inst);
void add_asymmetric_variables(CPXENVptr env, CPXLPptr lp, instance inst);

void add_symmetric_constraints(CPXENVptr env, CPXLPptr lp, instance inst);
void add_asymmetric_constraints(CPXENVptr env, CPXLPptr lp, instance inst);

void add_MTZ_sec(CPXENVptr env, CPXLPptr lp, instance inst);
void add_GG_sec(CPXENVptr env, CPXLPptr lp, instance inst);

double build_tsp_model(CPXENVptr env, CPXLPptr lp, instance inst, enum model_types model_type) {

	struct timeval start, end;
    gettimeofday(&start, NULL);

    switch (model_type) {
        case SYMMETRIC:
            add_symmetric_variables(env, lp, inst);
            add_symmetric_constraints(env, lp, inst);
            break;
        case SYMMETRIC_BENDERS:
            add_symmetric_variables(env, lp, inst);
            add_symmetric_constraints(env, lp, inst);
            break;
         case SYMMETRIC_BENDERS_CALLBACK:
            add_symmetric_variables(env, lp, inst);
            add_symmetric_constraints(env, lp, inst);
            break;

        case ASYMMETRIC_MTZ:
            add_asymmetric_variables(env, lp, inst);
            add_asymmetric_constraints(env, lp, inst);
            add_MTZ_sec(env, lp, inst);
            break;
        case ASYMMETRIC_GG:
            add_asymmetric_variables(env, lp, inst);
            add_asymmetric_constraints(env, lp, inst);
            add_GG_sec(env, lp, inst);
            break;

        default:
            print_error("unhandeled model type in variables");
    }

    /* save model file */
    char* filename;
    filename = (char*) calloc(100, sizeof(char));
    if (inst->model_folder == TSPLIB)    sprintf(filename, "../data_tsplib/%s/%s.", inst->model_name, inst->model_name);
    if (inst->model_folder == GENERATED) sprintf(filename, "../data_generated/%s/%s.", inst->model_name, inst->model_name);
    /* filename depend on model type */
    switch (model_type) {
        case SYMMETRIC:
            strcat(filename, "symmetric.lp");
            break;
        case SYMMETRIC_BENDERS:
            strcat(filename, "benders.lp");
            break;
        case ASYMMETRIC_MTZ:
            strcat(filename, "mtz.lp");
            break;
        case ASYMMETRIC_GG:
            strcat(filename, "gg.lp");

            break;
        default:
            break;
    }
    CPXwriteprob(env, lp, filename, NULL);
    free(filename);

    /* save build time */
    gettimeofday(&end, NULL);
    long seconds = (end.tv_sec - start.tv_sec);
    long micros = ((seconds * 1000000) + end.tv_usec) - (start.tv_usec);
	return micros / 1000.0;
}

void add_symmetric_variables(CPXENVptr env, CPXLPptr lp, instance inst) {
    /* varname ('C': continuous or 'B': binary)
     * lower and upper bound: lb <= x <= ub */
    char binary = 'B';
    double lb = 0.0;
    double ub = 1.0;

	/* ColumnNAME: array of array used to inject variables in CPLEX */
    char** cname = (char**) calloc(1, sizeof(char*));
    cname[0] = (char*) calloc(100, sizeof(char));

	/* add a single binary variables x(i,j) for i < j at the time */
    for (int i=0; i<inst->num_nodes; i++) for (int j=i+1; j<inst->num_nodes; j++) {

        /* write name of 1-indexed variable insede CPLEX
         * compute cost of var as distance x(i,j) */
        sprintf(cname[0], "x(%d,%d)", i+1, j+1);
        double obj = dist(i, j, inst);

		/* inject variable and test it's position (xpos) inside CPLEX */
        if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname)) {
            print_error("wrong CPXnewcols on x var.s");
        }
        if (CPXgetnumcols(env, lp)-1 != xpos(i, j, inst->num_nodes)) {
            print_error("wrong position for x var.s");
        }
    }
}

void add_asymmetric_variables(CPXENVptr env, CPXLPptr lp, instance inst) {
    /* varname ('C': continuous or 'B': binary)
     * upper bound 1.0: symmetric variable */
    char binary = 'B';
    double lb = 0.0;
    double ub = 1.0;

	/* ColumnNAME: array of array used to inject variables in CPLEX */
    char** cname = (char**) calloc(1, sizeof(char*));
    cname[0] = (char*) calloc(100, sizeof(char));

	/* add a single binary variables x(i,j) forall couples i,j at the time */
    for (int i=0; i<inst->num_nodes; i++) for (int j=0; j<inst->num_nodes; j++) {
        /* no self loops pls */
        if (i == j) ub = 0.0;
        else        ub = 1.0;

        /* write name of 1-indexed variable insede CPLEX
         * compute cost of var as distance x(i,j) */
        sprintf(cname[0], "x(%d,%d)", i+1, j+1);
        double obj = dist(i, j, inst);

		/* inject variable and test it's position (xxpos) inside CPLEX */
        if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname)) {
            print_error("wrong CPXnewcols on x var.s\n");
        }
        if (CPXgetnumcols(env, lp)-1 != xxpos(i, j, inst->num_nodes)) {
            print_error("wrong position for x var.s\n");
        }
    }
}

void add_symmetric_constraints(CPXENVptr env, CPXLPptr lp, instance inst) {
	/* right hand site scalar
	 * sense of the constraint: (E: equality, G: greater or equal, L: less or equal*/
    double rhs = 2.0;
    char sense = 'E';

	/* ColumnNAME: array of array used to inject variables in CPLEX */
    char** cname = (char**) calloc(1, sizeof(char*));
    cname[0] = (char*) calloc(100, sizeof(char));

	/* add one constraint at the time */
    for (int h=0; h<inst->num_nodes; h++) {
		/* fetch last row position in current model (better: the new one)
		 * write the constraint name inside CPLEX */
        int lastrow = CPXgetnumrows(env, lp);
        sprintf(cname[0], "degree(%d)", h+1);

		/* add the new constraint (with coefficent zero, null lhs) */
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
            print_error("wrong CPXnewrows [degree]");
        }
		/* and change the coefficent if node i is adiacent to h */
        for (int i=0; i<inst->num_nodes; i++) {
            if (i == h) continue;

            /* change the coefficent from 0 to 1.0 */
            if (CPXchgcoef(env, lp, lastrow, xpos(i, h, inst->num_nodes), 1.0)) {
                print_error("wrong CPXchgcoef [degree]");
            }
        }
    }
}

void add_asymmetric_constraints(CPXENVptr env, CPXLPptr lp, instance inst) {
	/* right hand site scalar
	 * sense of the constraint: (E: equality, G: greater or equal, L: less or equal*/
    double rhs = 1.0;
    char sense = 'E';

	/* ColumnNAME: array of array used to inject variables in CPLEX */
    char** cname = (char**) calloc(1, sizeof(char*));
    cname[0] = (char*) calloc(100, sizeof(char));

	/* add one constraint at the time */
    for (int h=0; h<inst->num_nodes; h++) {
		/* fetch last row position in current model (better: the new one)
		 * write the constraint name inside CPLEX */
        int lastrow = CPXgetnumrows(env, lp);
        sprintf(cname[0], "degree(%d)", h+1);

		/* add the new constraint (with coefficent zero, null lhs) */
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
            print_error("wrong CPXnewrows [degree]");
        }
		/* and change the coefficent if node i is adiacent to h */
        for (int i = 0; i<inst->num_nodes; i++) {
            if (i == h) continue;

            /* change the coefficent from 0 to 1.0 */
            if (CPXchgcoef(env, lp, lastrow, xxpos(i, h, inst->num_nodes), 1.0)) {
                print_error("wrong CPXchgcoef [degree]");
            }
        }
    }

    /* do the same but switch indexes (rhs = 1!) */
    for (int i=0; i<inst->num_nodes; i++) {
		/* fetch last row position in current model (better: the new one)
		 * write the constraint name inside CPLEX */
        int lastrow = CPXgetnumrows(env, lp);
        sprintf(cname[0], "degree(%d)", inst->num_nodes + i+1);

        /* add the new constraint (with coefficent zero, null lhs) */
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
            print_error("wrong CPXnewrows [degree]");
        }
        /* and change the coefficent if node i is adiacent to h */
        for (int h = 0; h<inst->num_nodes; h++) {
            /* the graph is complete so we only skip self loops */
            if (i == h) continue;

            /* change the coefficent from 0 to 1.0 */
            if (CPXchgcoef(env, lp, lastrow, xxpos(i, h, inst->num_nodes), 1.0)) {
                print_error("wrong CPXchgcoef [degree]");
            }
        }
    }
}

void add_MTZ_sec(CPXENVptr env, CPXLPptr lp, instance inst) {
	/* new continuous variables (integer does not matter here!)
	 * upper bound equal to m-1 for big-M constraint */
	char continuous = 'I';
	double lb = 1.0;
	double ub = inst->num_nodes-1;

	/* ColumnNAME: array of array used to inject variables in CPLEX */
	char** cname = (char**) calloc(1, sizeof(char*));
	cname[0] = (char*) calloc(100, sizeof(char));

    /* new variables associated with nodes, only n */
	for (int i=1; i<inst->num_nodes; i++) {

        /* cost is zero for new variables, they matter for new constraints only */
		sprintf(cname[0], "u(%d)", i+1);
		double obj = 0;

		/* inject variable and test it's position (upos) inside CPLEX */
		if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &continuous, cname)) {
			print_error("wrong CPXnewcols on x var.s");
		}
		if (CPXgetnumcols(env, lp)-1 != upos(i, inst->num_nodes)) {
			print_error(" wrong position for x var.s");
		}
	}


	int izero = 0;
	int index[3];
	double value[3];
	double big_M = inst->num_nodes - 1.0;
	double rhs = big_M -1.0;
	char sense = 'L';
	int nnz = 3;

	/* add lazy constraints  1.0 * u_i - 1.0 * u_j + M * x_ij <= M - 1
     * for each arc (i,j) not touching node 0 */
	for (int i=1; i<inst->num_nodes; i++) for (int j=1; j<inst->num_nodes; j++) {
        if (i == j) continue;

        sprintf(cname[0], "u-consistency for arc (%d,%d)", i+1, j+1);

        /* build constraint equation */
        index[0] = upos(i, inst->num_nodes);
        value[0] = 1.0;
        index[1] = upos(j, inst->num_nodes);
        value[1] = -1.0;
        index[2] = xxpos(i, j, inst->num_nodes);
        value[2] = big_M;

        if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index, value, cname) ) {
            print_error("wrong CPXlazyconstraints() for u-consistency");
        }
	}

	/* add lazy constraints 1.0 * x_ij + 1.0 * x_ji <= 1
     * for each arc (i,j) with i < j */
	rhs = 1.0;
	nnz = 2;
	for (int i=0; i<inst->num_nodes; i++) for (int j=i+1; j<inst->num_nodes; j++) {
        sprintf(cname[0], "SEC on node pair (%d,%d)", i+1, j+1);

        /* build constraint equation */
        index[0] = xxpos(i, j, inst->num_nodes);
        value[0] = 1.0;
        index[1] = xxpos(j, i, inst->num_nodes);
        value[1] = 1.0;
        if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index, value, cname)) {
            print_error("wrong CPXlazyconstraints on 2-node SECs");
        }
	}
}

void add_GG_sec(CPXENVptr env, CPXLPptr lp, instance inst) {
	/* new integer variables
	 * upper bound equal to m-1 */
	char integer = 'I';
	double lb = 0.0;
	double ub = inst->num_nodes-1;

	/* ColumnNAME: array of array used to inject variables in CPLEX */
	char** cname = (char**) calloc(1, sizeof(char*));
	cname[0] = (char*) calloc(100, sizeof(char));

	/* add a single flux variable y(i,j) forall i,j */
	for (int i=0; i<inst->num_nodes; i++) for(int j=0; j<inst->num_nodes; j++){
		if(i==j) continue;

        /* cost is zero for new variables, they matter for new constraints only */
		sprintf(cname[0], "y(%d)(%d)", i+1,j+1);
		double obj = 0;

		/* inject variable and test it's position (upos) inside CPLEX */
		if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &integer, cname)) {
			print_error("wrong CPXnewcols on x var.s");
		}
		if (CPXgetnumcols(env, lp)-1 != ypos(i, j, inst->num_nodes)) {
			print_error(" wrong position for x var.s");
		}
	}

	double rhs = 0.0;
	char sense = 'L';

	for (int i=0; i<inst->num_nodes; i++) for (int h=0; h<inst->num_nodes; h++) {
        if (i == h) continue;

        /* fetch last row position in current model (better: the new one)
         * in CPLEX row == constraint */
        int lastrow = CPXgetnumrows(env, lp);
        sprintf(cname[0], "xy(%d)(%d)", i+1,h+1);

        /* add the new constraint and change coeff from zero to 1 and 1-m */
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
            print_error("wrong CPXnewrows [degree]");
        }
        if (CPXchgcoef(env, lp, lastrow, xxpos(i, h, inst->num_nodes), 1-inst->num_nodes)) {
            print_error("wrong CPXchgcoef [degree]");
        }
        if (CPXchgcoef(env, lp, lastrow, ypos(i, h, inst->num_nodes), 1.0)) {
            print_error("wrong CPXchgcoef [degree]");
        }
    }

    rhs = inst->num_nodes - 1;
    sense = 'E';
    int lastrow = CPXgetnumrows(env, lp);
    sprintf(cname[0], "y1j");

    /* add the new constraint (with coefficent zero, null lhs) */
    if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
        print_error("wrong CPXnewrows [degree]");
    }

    for (int h=1; h<inst->num_nodes; h++) {
        if (CPXchgcoef(env, lp, lastrow, ypos(0, h, inst->num_nodes), 1.0)) {
            print_error("wrong CPXchgcoef [degree]");
        }
    }

    rhs = 1;
    sense = 'E';
    for (int j=1; j<inst->num_nodes; j++) {

        int lastrow = CPXgetnumrows(env, lp);
        sprintf(cname[0], "y(%d)", j+1);

        /* add the new constraint (with coefficent zero, null lhs) */
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
            print_error("wrong CPXnewrows [degree]");
        }

        /* and change the coefficent if node i is adiacent to h */
        for (int i = 0; i<inst->num_nodes; i++) {
            /* the graph is complete so we only skip self loops */
            if (i == j) continue;

            /* change the coefficent from 0 to 1.0 */
            if (CPXchgcoef(env, lp, lastrow, ypos(i, j, inst->num_nodes), 1.0)) {
                print_error("wrong CPXchgcoef [degree]");
            }
        }

        /* and change the coefficent if node i is adiacent to h */
        for (int h = 0; h<inst->num_nodes; h++) {
            /* the graph is complete so we only skip self loops */
            if (h == j) continue;

            /* change the coefficent from 0 to 1.0 */
            if (CPXchgcoef(env, lp, lastrow, ypos(j, h, inst->num_nodes), -1.0)) {
                print_error("wrong CPXchgcoef [degree]");
            }
        }
    }
}

void add_BENDERS_sec(CPXENVptr env, CPXLPptr lp, solution sol) {
    int* visited = (int*) calloc(sol->num_edges, sizeof(int));

    /* ColumnNAME: array of array used to inject variables in CPLEX */
    char** cname = (char**) calloc(1, sizeof(char*));
    cname[0] = (char*) calloc(100, sizeof(char));

    for (int i=0; i<sol->num_edges; i++) {
        /* i could be part of a subtour already visited */
        if (visited[i]) continue;

        double rhs = 0.0;
        char sense = 'L';

        /* compute the rhs: number of nodes in a subtour, unknown a priori */
        int j = i;
        while ((j = sol->parent[j]) != i) rhs += 1.0;


		/* fetch last row position in current model (better: the new one)
		 * write the constraint name inside CPLEX */
        int lastrow = CPXgetnumrows(env, lp);
        sprintf(cname[0], "benders_sec(%d)", i+1);

        /* add the new constraint (with coefficent zero, null lhs) */
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
            print_error("wrong CPXnewrows [degree]");
        }

        /* i: iterate over subtours, start element
         * j: iterate inside subtour, first index
         * k: iterate inside subtour, second index
         * -> we will add x_jk to the subtour
         */
        j = i;
        do {
            int k = sol->parent[j];
            do {
                /* change the coefficent from 0 to 1.0 */
                if (CPXchgcoef(env, lp, lastrow, xpos(j, k, sol->num_edges), 1.0)) {
                    print_error("wrong CPXchgcoef [degree]");
                }
                k = sol->parent[k];
            } while (k != i);

            j = sol->parent[j];
        } while (sol->parent[j] != i);

        /* visit all the nodes in the set */
        j = i;
        while ((j = sol->parent[j]) != i) visited[j] = 1;
    }

    free(visited);
}

int CPXPUBLIC add_BENDERS_sec_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle) {

    solution sol = (solution) userhandle;

    /* number of columns is n chooses 2 */
    int num_cols = sol->num_edges * (sol->num_edges-1) /2;

    double* xstar = (double*) calloc(num_cols, sizeof(double));
    double objval = CPX_INFBOUND;

    if (CPXcallbackgetcandidatepoint(context, xstar, 0, num_cols - 1, &objval)) {
        print_error("CPXcallbackgetcandidatepoint error");
    }

    /* retreive sol and fill sol's internal parent vector */
    int* visited = (int*) calloc(sol->num_edges, sizeof(int));
    retreive_symmetric_solution(xstar, sol);


    for (int i=0; i<sol->num_edges; i++) {
        /* i could be part of a subtour already visited */
        if (visited[i]) continue;

        int nnz = 0;
        int izero = 0;
        double rhs = 0.0;
        char sense = 'L';
        int* index = (int*) calloc(num_cols, sizeof(int));
        double* value = (double*) calloc(num_cols, sizeof(double));


        /* compute the rhs: number of nodes in a subtour, unknown a priori */
        int j = i;
        while ((j = sol->parent[j]) != i) rhs += 1.0;

        /* i: iterate over subtours, start element
         * j: iterate inside subtour, first index
         * k: iterate inside subtour, second index
         * -> we will add x_jk to the subtour
         */
        j = i;
        do {
            int k = sol->parent[j];
            do {
                index[nnz] = xpos(j, k, sol->num_edges);
                value[nnz++] = 1.0;
                k = sol->parent[k];
            } while (k != i);

            j = sol->parent[j];
        } while (sol->parent[j] != i);

        /* visit all the nodes in the set */
        j = i;
        while ((j = sol->parent[j]) != i) visited[j] = 1;

        /* finally set the callback for rejecte the incumbent */
        if (rhs != sol->num_edges - 1) {
        	if (CPXcallbackrejectcandidate(context, 1, nnz, &rhs, &sense, &izero, index, value)) {
            	print_error("CPXcallbackrejectcandidate() error");
            }
        }

        free(index);
        free(value);
    }

    free(visited);

    return 0;
}
