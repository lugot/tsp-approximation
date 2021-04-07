#include "tsp.h"
#include "utils.h"
#include "globals.h"
#include "string.h"
#include <stdio.h>
#include <stdlib.h>
#include <cplex.h>
#include <sys/time.h>

void add_symmetric_variables(instance inst, CPXENVptr env, CPXLPptr lp);
void add_asymmetric_variables(instance inst, CPXENVptr env, CPXLPptr lp);

void add_symmetric_constraints(instance inst, CPXENVptr env, CPXLPptr lp);
void add_asymmetric_constraints(instance inst, CPXENVptr env, CPXLPptr lp);

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

	struct timeval start, end;
    gettimeofday(&start, NULL);

    switch (model_type) {
        case SYMMETRIC:
            add_symmetric_variables(inst, env, lp);
            add_symmetric_constraints(inst, env, lp);
            break;
        case SYMMETRIC_BENDERS:
            add_symmetric_variables(inst, env, lp);
            add_symmetric_constraints(inst, env, lp);
            break;

        case ASYMMETRIC_MTZ:
            add_asymmetric_variables(inst, env, lp);
            add_asymmetric_constraints(inst, env, lp);
            add_MTZ_subtour(inst, env, lp);
            break;
        case ASYMMETRIC_GG:
            add_asymmetric_variables(inst, env, lp);
            add_asymmetric_constraints(inst, env, lp);
            add_GG_subtour(inst, env, lp);
            break;

        default:
            print_error("unhandeled model type in variables");
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

void add_symmetric_variables(instance inst, CPXENVptr env, CPXLPptr lp) {
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
        if (CPXgetnumcols(env, lp)-1 != xpos(i, j, inst)) {
            print_error("wrong position for x var.s");
        }
    }

    free(cname[0]);
    free(cname);
}

void add_asymmetric_variables(instance inst, CPXENVptr env, CPXLPptr lp) {
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
        if (CPXgetnumcols(env, lp)-1 != xxpos(i, j, inst)) {
            print_error("wrong position for x var.s\n");
        }
    }

    free(cname[0]);
    free(cname);
}

void add_symmetric_constraints(instance inst, CPXENVptr env, CPXLPptr lp) {
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
            if (CPXchgcoef(env, lp, lastrow, xpos(i, h, inst), 1.0)) {
                print_error("wrong CPXchgcoef [degree]");
            }
        }
    }

    free(cname[0]);
    free(cname);
}

void add_asymmetric_constraints(instance inst, CPXENVptr env, CPXLPptr lp) {
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
            if (CPXchgcoef(env, lp, lastrow, xxpos(i, h, inst), 1.0)) {
                print_error("wrong CPXchgcoef [degree]");
            }
        }
    }

    /* do the same but switch indexes (rhs = 1!) */
    for (size_t i=0; i<inst->num_nodes; i++) {
		/* fetch last row position in current model (better: the new one)
		 * write the constraint name inside CPLEX */
        int lastrow = CPXgetnumrows(env, lp);
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

    free(cname[0]);
    free(cname);
}

void add_MTZ_subtour(instance inst, CPXENVptr env, CPXLPptr lp) {
	/* new continuous variables (integer does not matter here!)
	 * upper bound equal to m-1 for big-M constraint */
	char continuous = 'I';
	double lb = 1.0;
	double ub = inst->num_nodes-1;

	/* ColumnNAME: array of array used to inject variables in CPLEX */
	char** cname = (char**) calloc(1, sizeof(char*));
	cname[0] = (char*) calloc(100, sizeof(char));

    /* new variables associated with nodes, only n */
	for (size_t i=1; i<inst->num_nodes; i++) {

        /* cost is zero for new variables, they matter for new constraints only */
		sprintf(cname[0], "u(%ld)", i+1);
		double obj = 0;

		/* inject variable and test it's position (upos) inside CPLEX */
		if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &continuous, cname)) {
			print_error("wrong CPXnewcols on x var.s");
		}
		if (CPXgetnumcols(env, lp)-1 != upos(i, inst)) {
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
        index[0] = upos(i,inst);
        value[0] = 1.0;
        index[1] = upos(j,inst);
        value[1] = -1.0;
        index[2] = xxpos(i,j,inst);
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
        index[0] = xxpos(i,j,inst);
        value[0] = 1.0;
        index[1] = xxpos(j,i,inst);
        value[1] = 1.0;
        if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index, value, cname)) {
            print_error("wrong CPXlazyconstraints on 2-node SECs");
        }
	}

    free(cname[0]);
    free(cname);
}

void add_GG_subtour(instance inst, CPXENVptr env, CPXLPptr lp) {
	/* new integer variables
	 * upper bound equal to m-1 */
	char integer = 'I';
	double lb = 0.0;
	double ub = inst->num_nodes-1;

	/* ColumnNAME: array of array used to inject variables in CPLEX */
	char** cname = (char**) calloc(1, sizeof(char*));
	cname[0] = (char*) calloc(100, sizeof(char));

	/* add a single flux variable y(i,j) forall i,j */
	for (size_t i=0; i<inst->num_nodes; i++) for(size_t j=0; j<inst->num_nodes; j++){
		if(i==j) continue;

        /* cost is zero for new variables, they matter for new constraints only */
		sprintf(cname[0], "y(%ld)(%ld)", i+1,j+1);
		double obj = 0;

		/* inject variable and test it's position (upos) inside CPLEX */
		if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &integer, cname)) {
			print_error("wrong CPXnewcols on x var.s");
		}
		if (CPXgetnumcols(env, lp)-1 != ypos(i, j, inst)) {
			print_error(" wrong position for x var.s");
		}
	}

	double rhs = 0;
	char sense = 'L';

	for (size_t i=0; i<inst->num_nodes; i++) for (size_t h=0; h<inst->num_nodes; h++) {
        if (i == h) continue;

        /* fetch last row position in current model (better: the new one)
         * in CPLEX row == constraint */
        int lastrow = CPXgetnumrows(env, lp);
        sprintf(cname[0], "xy(%ld)(%ld)", i+1,h+1);

        /* add the new constraint and change coeff from zero to 1 and 1-m */
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
            print_error("wrong CPXnewrows [degree]");
        }
        if (CPXchgcoef(env, lp, lastrow, xxpos(i, h, inst), 1-inst->num_nodes)) {
            print_error("wrong CPXchgcoef [degree]");
        }
        if (CPXchgcoef(env, lp, lastrow, ypos(i, h, inst), 1.0)) {
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

    for (size_t h=1; h<inst->num_nodes; h++) {
        if (CPXchgcoef(env, lp, lastrow, ypos(0, h, inst), 1.0)) {
            print_error("wrong CPXchgcoef [degree]");
        }
    }

    rhs = 1;
    sense = 'E';
    for (size_t j=1; j<inst->num_nodes; j++) {

        int lastrow = CPXgetnumrows(env, lp);
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

    free(cname[0]);
    free(cname);
}
<<<<<<< HEAD

void add_cool_subtour_elimination(instance inst, CPXENVptr env, CPXLPptr lp) {
    int visited[inst->num_nodes];
    for (int i=0; i<inst->num_nodes; i++) visited[i]=0;
    /*memset(visited, 0, inst->num_nodes*sizeof(int));*/

    /* ColumnNAME: array of array used to inject variables in CPLEX */
    char** cname = (char**) calloc(1, sizeof(char*));
    cname[0] = (char*) calloc(100, sizeof(char));

    for (int i=0; i<inst->num_nodes; i++) {
        if (visited[i]) continue;

        int act_set = uf_find_set(uf, i);
        int set_size = uf_set_size(uf, i);

        if (set_size <= 2) continue;

        double rhs = uf_set_size(uf, i) - 1;
        char sense = 'L';

		/* fetch last row position in current model (better: the new one)
		 * write the constraint name inside CPLEX */
        int lastrow = CPXgetnumrows(env, lp);
        sprintf(cname[0], "subelimination(%d)", i+1);

        /* add the new constraint (with coefficent zero, null lhs) */
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
            print_error("wrong CPXnewrows [degree]");
        }

        // TODO: this is dumb af and overcomplicated

        int j=i;
        do {
            int k = uf->next[j];
            do {
                /* change the coefficent from 0 to 1.0 */
                if (CPXchgcoef(env, lp, lastrow, xpos(j, k, inst), 1.0)) {
                    print_error("wrong CPXchgcoef [degree]");
                }
                k = uf->next[k];
            } while (k != i);

            j = uf->next[j];
        } while (uf->next[j] != i);

        /* visit all the nodes in the set */
        for (int j=0; j<inst->num_nodes; j++) {
            if (uf_find_set(uf, j) == act_set) visited[j] = 1;
        }
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
=======
>>>>>>> parent of f2d3960 (added benders)
