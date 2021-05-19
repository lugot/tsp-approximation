#include "../include/model_builder.h"

#include <assert.h>
#include <concorde.h>
#include <cplex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

#include "../include/adjlist.h"
#include "../include/globals.h"
#include "../include/solvers.h"
#include "../include/tsp.h"
#include "../include/union_find.h"
#include "../include/utils.h"

void add_symm_variables(CPXENVptr env, CPXLPptr lp, instance inst);
void add_asymm_variables(CPXENVptr env, CPXLPptr lp, instance inst);

void add_symm_constraints(CPXENVptr env, CPXLPptr lp, instance inst);
void add_asymm_constraints(CPXENVptr env, CPXLPptr lp, instance inst);

void add_MTZ_static_sec(CPXENVptr env, CPXLPptr lp, instance inst);
void add_MTZ_lazy_sec(CPXENVptr env, CPXLPptr lp, instance inst);
void add_GGlit_static_sec(CPXENVptr env, CPXLPptr lp, instance inst);
void add_GGlect_static_sec(CPXENVptr env, CPXLPptr lp, instance inst);
void add_GGlit_lazy_sec(CPXENVptr env, CPXLPptr lp, instance inst);

int CPXPUBLIC add_BENDERS_sec_callback_candidate(CPXCALLBACKCONTEXTptr context,
                                                 solution sol);
int CPXPUBLIC add_BENDERS_sec_callback_relaxation(CPXCALLBACKCONTEXTptr context,
                                                  solution sol);
int doit_fn_concorde(double cutval, int cutcount, int* cut, void* in);

double build_tsp_model(CPXENVptr env, CPXLPptr lp, instance inst,
                       enum model_types model_type) {
    struct timeval start, end;
    gettimeofday(&start, NULL);

    switch (model_type) {
        case NOSEC:
            add_symm_variables(env, lp, inst);
            add_symm_constraints(env, lp, inst);
            break;
        case BENDERS:
            add_symm_variables(env, lp, inst);
            add_symm_constraints(env, lp, inst);
            break;
        case BENDERS_CALLBACK:
            add_symm_variables(env, lp, inst);
            add_symm_constraints(env, lp, inst);
            break;
        case HARD_FIXING:
            add_symm_variables(env, lp, inst);
            add_symm_constraints(env, lp, inst);
            break;
        case SOFT_FIXING:
            add_symm_variables(env, lp, inst);
            add_symm_constraints(env, lp, inst);
            break;
        case MTZ_STATIC:
            add_asymm_variables(env, lp, inst);
            add_asymm_constraints(env, lp, inst);
            add_MTZ_static_sec(env, lp, inst);
            break;
        case MTZ_LAZY:
            add_asymm_variables(env, lp, inst);
            add_asymm_constraints(env, lp, inst);
            add_MTZ_lazy_sec(env, lp, inst);
            break;
        case GGLIT_STATIC:
            add_asymm_variables(env, lp, inst);
            add_asymm_constraints(env, lp, inst);
            add_GGlit_static_sec(env, lp, inst);
            break;
        case GGLECT_STATIC:
            add_asymm_variables(env, lp, inst);
            add_asymm_constraints(env, lp, inst);
            add_GGlect_static_sec(env, lp, inst);
            break;
        case GGLIT_LAZY:
            add_asymm_variables(env, lp, inst);
            add_asymm_constraints(env, lp, inst);
            add_GGlit_lazy_sec(env, lp, inst);
            break;

        case MST:
        case GRASP:
        case GREEDY:
        case EXTRA_MILEAGE:
        case OPTIMAL_TOUR:
        case VNS:
            print_error("unhandeled model type in variables");
    }

    /* save model file */
    char* filename;
    int bufsize = 100;
    filename = (char*)calloc(bufsize, sizeof(char));
    if (inst->model_folder == TSPLIB)
        snprintf(filename, bufsize, "../data_tsplib/%s/%s.", inst->model_name,
                 inst->model_name);
    if (inst->model_folder == GENERATED)
        snprintf(filename, bufsize, "../data_generated/%s/%s.",
                 inst->model_name, inst->model_name);

    /* filename depend on model type */
    char* model_type_str = model_type_tostring(model_type);
    snprintf(filename + strlen(filename), bufsize, "%s.lp", model_type_str);
    CPXwriteprob(env, lp, filename, "LP");

    /* save build time */
    gettimeofday(&end, NULL);
    int64_t seconds = (end.tv_sec - start.tv_sec);
    int64_t micros = ((seconds * 1000000) + end.tv_usec) - (start.tv_usec);
    return micros / 1000.0;

    free(model_type_str);
    free(filename);
}

void add_symm_variables(CPXENVptr env, CPXLPptr lp, instance inst) {
    /* varname ('C': continuous or 'B': binary)
     * lower and upper bound: lb <= x <= ub */
    char binary = 'B';
    double lb = 0.0;
    double ub = 1.0;

    /* ColumnNAME: array of array used to inject variables in CPLEX */
    char** cname = (char**)calloc(1, sizeof(char*));
    cname[0] = (char*)calloc(100, sizeof(char));

    /* add a single binary variables x(i,j) for i < j at the time */
    int nnodes = inst->nnodes;
    for (int i = 0; i < nnodes; i++) {
        for (int j = i + 1; j < nnodes; j++) {
            /* write name of 1-indexed variable insede CPLEX
             * compute cost of var as distance x(i,j) */
            snprintf(cname[0], strlen(cname[0]), "x(%d-%d)", i + 1, j + 1);
            double obj = dist(i, j, inst);

            /* inject variable and test it's position (xpos) inside CPLEX */
            if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname)) {
                print_error("wrong CPXnewcols on x var.s");
            }
            if (CPXgetnumcols(env, lp) - 1 != xpos(i, j, nnodes)) {
                print_error("wrong position for x var.s");
            }
        }
    }

    /* add the number of CPLEX columns to the instance */
    /* int this case is n chooses 2 cause of symmetry */
    inst->ncols = nnodes * (nnodes - 1) / 2;

    free(cname[0]);
    free(cname);
}

void add_asymm_variables(CPXENVptr env, CPXLPptr lp, instance inst) {
    /* varname ('C': continuous or 'B': binary)
     * upper bound 1.0: symmetric variable */
    char binary = 'B';
    double lb = 0.0;
    double ub = 1.0;

    /* ColumnNAME: array of array used to inject variables in CPLEX */
    char** cname = (char**)calloc(1, sizeof(char*));
    cname[0] = (char*)calloc(100, sizeof(char));

    /* add a single binary variables x(i,j) forall couples i,j at the time */
    int nnodes = inst->nnodes;
    for (int i = 0; i < nnodes; i++) {
        for (int j = 0; j < nnodes; j++) {
            /* no self loops pls */
            if (i == j)
                ub = 0.0;
            else
                ub = 1.0;

            /* write name of 1-indexed variable insede CPLEX
             * compute cost of var as distance x(i,j) */
            snprintf(cname[0], strlen(cname[0]), "x(%d-%d)", i + 1, j + 1);
            double obj = dist(i, j, inst);

            /* inject variable and test it's position (xxpos) inside CPLEX */
            if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname)) {
                print_error("wrong CPXnewcols on x var.s\n");
            }
            if (CPXgetnumcols(env, lp) - 1 != xxpos(i, j, nnodes)) {
                print_error("wrong position for x var.s\n");
            }
        }
    }

    /* add the number of CPLEX columns to the instance */
    /* int this case is n * (n - 1) cause of asymmetry */
    inst->ncols = nnodes * (nnodes - 1);

    free(cname[0]);
    free(cname);
}

void add_symm_constraints(CPXENVptr env, CPXLPptr lp, instance inst) {
    /* right hand site scalar
     * sense of the constraint: (E: equality, G: greater or equal, L: less or
     * equal*/
    double rhs = 2.0;
    char sense = 'E';

    /* ColumnNAME: array of array used to inject variables in CPLEX */
    char** cname = (char**)calloc(1, sizeof(char*));
    cname[0] = (char*)calloc(100, sizeof(char));

    /* add one constraint at the time */ /*{{{*/ /*}}}*/
    int nnodes = inst->nnodes;
    for (int h = 0; h < nnodes; h++) {
        /* fetch last row position in current model (better: the new one)
         * write the constraint name inside CPLEX */
        int lastrow = CPXgetnumrows(env, lp);
        snprintf(cname[0], strlen(cname[0]), "degree(%d)", h + 1);

        /* add the new constraint (with coefficent zero, null lhs) */
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
            print_error("wrong CPXnewrows [degree]");
        }
        /* and change the coefficent if node i is adiacent to h */
        for (int i = 0; i < nnodes; i++) {
            if (i == h) continue;

            /* change the coefficent from 0 to 1.0 */
            if (CPXchgcoef(env, lp, lastrow, xpos(i, h, nnodes), 1.0)) {
                print_error("wrong CPXchgcoef [degree]");
            }
        }
    }

    free(cname[0]);
    free(cname);
}

void add_asymm_constraints(CPXENVptr env, CPXLPptr lp, instance inst) {
    /* right hand site scalar
     * sense of the constraint: (E: equality, G: greater or equal, L: less or
     * equal*/
    double rhs = 1.0;
    char sense = 'E';

    /* ColumnNAME: array of array used to inject variables in CPLEX */
    char** cname = (char**)calloc(1, sizeof(char*));
    cname[0] = (char*)calloc(100, sizeof(char));

    /* add one constraint at the time */
    int nnodes = inst->nnodes;
    for (int h = 0; h < nnodes; h++) {
        /* fetch last row position in current model (better: the new one)
         * write the constraint name inside CPLEX */
        int lastrow = CPXgetnumrows(env, lp);
        snprintf(cname[0], strlen(cname[0]), "degree(%d)", h + 1);

        /* add the new constraint (with coefficent zero, null lhs) */
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
            print_error("wrong CPXnewrows [degree]");
        }
        /* and change the coefficent if node i is adiacent to h */
        for (int i = 0; i < nnodes; i++) {
            if (i == h) continue;

            /* change the coefficent from 0 to 1.0 */
            if (CPXchgcoef(env, lp, lastrow, xxpos(i, h, nnodes), 1.0)) {
                print_error("wrong CPXchgcoef [degree]");
            }
        }
    }

    /* do the same but switch indexes (rhs = 1!) */
    for (int i = 0; i < nnodes; i++) {
        /* fetch last row position in current model (better: the new one)
         * write the constraint name inside CPLEX */
        int lastrow = CPXgetnumrows(env, lp);
        snprintf(cname[0], strlen(cname[0]), "degree(%d)", nnodes + i + 1);

        /* add the new constraint (with coefficent zero, null lhs) */
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
            print_error("wrong CPXnewrows [degree]");
        }
        /* and change the coefficent if node i is adiacent to h */
        for (int h = 0; h < nnodes; h++) {
            /* the graph is complete so we only skip self loops */
            if (i == h) continue;

            /* change the coefficent from 0 to 1.0 */
            if (CPXchgcoef(env, lp, lastrow, xxpos(i, h, nnodes), 1.0)) {
                print_error("wrong CPXchgcoef [degree]");
            }
        }
    }

    free(cname[0]);
    free(cname);
}

void add_MTZ_static_sec(CPXENVptr env, CPXLPptr lp, instance inst) {
    int nnodes = inst->nnodes;
    /* new continuous variables (integer does not matter here!)
     * upper bound equal to m-1 for big-M constraint */
    char continuous = 'I';
    double lb = 1.0;
    double ub = nnodes - 1;

    /* ColumnNAME: array of array used to inject variables in CPLEX */
    char** cname = (char**)calloc(1, sizeof(char*));
    cname[0] = (char*)calloc(100, sizeof(char));

    /* new variables associated with nodes, only n */
    for (int i = 1; i < nnodes; i++) {
        /* cost is zero for new variables, they matter for new constraints only
         */
        snprintf(cname[0], strlen(cname[0]), "u(%d)", i + 1);
        double obj = 0;

        /* inject variable and test it's position (upos) inside CPLEX */
        if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &continuous, cname)) {
            print_error("wrong CPXnewcols on x var.s");
        }
        if (CPXgetnumcols(env, lp) - 1 != upos(i, nnodes)) {
            print_error(" wrong position for x var.s");
        }
    }

    /* update the number of columns */
    inst->ncols += nnodes;

    int izero = 0;
    int index[3];
    double value[3];
    double big_M = nnodes - 1.0;
    double rhs = big_M - 1.0;
    char sense = 'L';
    int nnz = 3;

    /* add static constraints  1.0 * u_i - 1.0 * u_j + M * x_ij <= M - 1
     * for each arc (i,j) not touching node 0 */
    for (int i = 1; i < nnodes; i++)
        for (int j = 1; j < nnodes; j++) {
            if (i == j) continue;
            /* fetch last row position in current model (better: the new one)
             * write the constraint name inside CPLEX */
            snprintf(cname[0], strlen(cname[0]),
                     "u-consistency for arc (%d,%d)", i + 1, j + 1);

            /* build constraint equation */
            index[0] = upos(i, nnodes);
            value[0] = 1.0;
            index[1] = upos(j, nnodes);
            value[1] = -1.0;
            index[2] = xxpos(i, j, nnodes);
            value[2] = big_M;

            if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index,
                           value, NULL, cname)) {
                print_error("wrong CPXlazyconstraints() for u-consistency");
            }
        }

    /* add static constraints 1.0 * x_ij + 1.0 * x_ji <= 1
     * for each arc (i,j) with i < j */
    rhs = 1.0;
    nnz = 2;
    for (int i = 0; i < nnodes; i++)
        for (int j = i + 1; j < nnodes; j++) {
            snprintf(cname[0], strlen(cname[0]), "SEC on node pair (%d,%d)",
                     i + 1, j + 1);

            /* build constraint equation */
            index[0] = xxpos(i, j, nnodes);
            value[0] = 1.0;
            index[1] = xxpos(j, i, nnodes);
            value[1] = 1.0;
            if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index,
                           value, NULL, cname)) {
                print_error("wrong CPXlazyconstraints on 2-node SECs");
            }
        }

    free(cname[0]);
    free(cname);
}

void add_MTZ_lazy_sec(CPXENVptr env, CPXLPptr lp, instance inst) {
    int nnodes = inst->nnodes;
    /* new continuous variables (integer does not matter here!)
     * upper bound equal to m-1 for big-M constraint */
    char continuous = 'I';
    double lb = 1.0;
    double ub = nnodes - 1;

    /* ColumnNAME: array of array used to inject variables in CPLEX */
    char** cname = (char**)calloc(1, sizeof(char*));
    cname[0] = (char*)calloc(100, sizeof(char));

    /* new variables associated with nodes, only n */
    for (int i = 1; i < nnodes; i++) {
        /* cost is zero for new variables, they matter for new constraints only
         */
        snprintf(cname[0], strlen(cname[0]), "u(%d)", i + 1);
        double obj = 0;

        /* inject variable and test it's position (upos) inside CPLEX */
        if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &continuous, cname)) {
            print_error("wrong CPXnewcols on x var.s");
        }
        if (CPXgetnumcols(env, lp) - 1 != upos(i, nnodes)) {
            print_error(" wrong position for x var.s");
        }
    }

    /* update the number of columns */
    inst->ncols += nnodes;

    int izero = 0;
    int index[3];
    double value[3];
    double big_M = nnodes - 1.0;
    double rhs = big_M - 1.0;
    char sense = 'L';
    int nnz = 3;

    /* add lazy constraints  1.0 * u_i - 1.0 * u_j + M * x_ij <= M - 1
     * for each arc (i,j) not touching node 0 */
    for (int i = 1; i < nnodes; i++) {
        for (int j = 1; j < nnodes; j++) {
            if (i == j) continue;

            snprintf(cname[0], strlen(cname[0]),
                     "u-consistency for arc (%d,%d)", i + 1, j + 1);

            /* build constraint equation */
            index[0] = upos(i, nnodes);
            value[0] = 1.0;
            index[1] = upos(j, nnodes);
            value[1] = -1.0;
            index[2] = xxpos(i, j, nnodes);
            value[2] = big_M;

            if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero,
                                      index, value, cname)) {
                print_error("wrong CPXlazyconstraints() for u-consistency");
            }
        }
    }

    /* add lazy constraints 1.0 * x_ij + 1.0 * x_ji <= 1
     * for each arc (i,j) with i < j */
    rhs = 1.0;
    nnz = 2;
    for (int i = 0; i < nnodes; i++)
        for (int j = i + 1; j < nnodes; j++) {
            snprintf(cname[0], strlen(cname[0]), "SEC on node pair (%d,%d)",
                     i + 1, j + 1);

            /* build constraint equation */
            index[0] = xxpos(i, j, nnodes);
            value[0] = 1.0;
            index[1] = xxpos(j, i, nnodes);
            value[1] = 1.0;
            if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero,
                                      index, value, cname)) {
                print_error("wrong CPXlazyconstraints on 2-node SECs");
            }
        }

    free(cname[0]);
    free(cname);
}

void add_GGlit_static_sec(CPXENVptr env, CPXLPptr lp, instance inst) {
    int nnodes = inst->nnodes;
    /* new integer variables
     * upper bound equal to m-1 */
    char integer = 'I';
    double lb = 0.0;
    double ub = nnodes - 1;

    /* ColumnNAME: array of array used to inject variables in CPLEX */
    char** cname = (char**)malloc(sizeof(char*));
    int bufsize = 100;
    cname[0] = (char*)calloc(bufsize, sizeof(char));

    /* add a single flux variable y(i,j) forall i,j */
    for (int i = 0; i < nnodes; i++) {
        for (int j = 0; j < nnodes; j++) {
            if (i == j) continue;

            /* cost is zero for new variables, they matter for new constraints
             * only */
            snprintf(cname[0], bufsize, "y(%d)(%d)", i + 1, j + 1);
            double obj = 0;

            /* inject variable and test it's position (upos) inside CPLEX */
            if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &integer, cname)) {
                print_error("wrong CPXnewcols on x var.s");
            }
            if (CPXgetnumcols(env, lp) - 1 != ypos(i, j, nnodes)) {
                print_error(" wrong position for x var.s");
            }
        }
    }

    /* update the number of columns */
    inst->ncols += nnodes * nnodes;

    double rhs = 0.0;
    char sense = 'L';
    /* linking constraint: y_ij <= (n-1) x_ij forall i,j in V
     * transformed in (1-n) x_ij + y_ij <= 0 */
    for (int i = 0; i < nnodes; i++) {
        for (int j = 0; j < nnodes; j++) {
            if (i == j) continue;

            /* tighten constraint */
            int x_coeff;
            if (i != 0)
                x_coeff = 2 - nnodes;
            else
                x_coeff = 1 - nnodes;

            /* fetch last row position in current model (better: the new one)
             * in CPLEX row == constraint */
            int lastrow = CPXgetnumrows(env, lp);
            snprintf(cname[0], bufsize, "xy(%d)(%d)", i + 1, j + 1);

            /* add the new constraint and change coeff from zero to 1 and 1-m */
            if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
                print_error("wrong CPXnewrows [degree]");
            }
            if (CPXchgcoef(env, lp, lastrow, xxpos(i, j, nnodes), x_coeff)) {
                print_error("wrong CPXchgcoef [degree]");
            }
            if (CPXchgcoef(env, lp, lastrow, ypos(i, j, nnodes), 1.0)) {
                print_error("wrong CPXchgcoef [degree]");
            }
        }
    }

    /* sum j in V \ {1} y_1j = n - 1 */
    rhs = nnodes - 1;
    sense = 'E';
    int lastrow = CPXgetnumrows(env, lp);
    snprintf(cname[0], bufsize, "y1j");

    /* add the new constraint (with coefficent zero, null lhs) */
    if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
        print_error("wrong CPXnewrows [degree]");
    }
    for (int j = 1; j < nnodes; j++) {
        if (CPXchgcoef(env, lp, lastrow, ypos(0, j, nnodes), 1.0)) {
            print_error("wrong CPXchgcoef [degree]");
        }
    }

    /* flow control:
     * sum i != j  y_ij - sum k != j y_jk = 1 forall j in V \ {1} */
    rhs = 1.0;
    sense = 'E';
    for (int j = 1; j < nnodes; j++) {
        lastrow = CPXgetnumrows(env, lp);
        snprintf(cname[0], bufsize, "y(%d)", j + 1);

        /* add the new constraint (with coefficent zero, null lhs) */
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
            print_error("wrong CPXnewrows [degree]");
        }

        /* and change the coefficent if node i is adiacent to j */
        for (int i = 0; i < nnodes; i++) {
            /* the graph is complete so we only skip self loops */
            if (i == j) continue;

            /* change the coefficent from 0 to 1.0 */
            if (CPXchgcoef(env, lp, lastrow, ypos(i, j, nnodes), 1.0)) {
                print_error("wrong CPXchgcoef [degree]");
            }
        }

        /* and change the coefficent if node j is adiacent to h */
        for (int k = 0; k < nnodes; k++) {
            /* the graph is complete so we only skip self loops */
            if (k == j) continue;

            /* change the coefficent from 0 to 1.0 */
            if (CPXchgcoef(env, lp, lastrow, ypos(j, k, nnodes), -1.0)) {
                print_error("wrong CPXchgcoef [degree]");
            }
        }
    }

    free(cname[0]);
    free(cname);
}

void add_GGlect_static_sec(CPXENVptr env, CPXLPptr lp, instance inst) {
    int nnodes = inst->nnodes;
    /* new integer variables
     * upper bound equal to m-1 */
    char integer = 'I';
    double lb = 0.0;
    double ub = nnodes - 2; /* usually this but set later */

    /* ColumnNAME: array of array used to inject variables in CPLEX */
    char** cname = (char**)calloc(1, sizeof(char*));
    int bufsize = 100;
    cname[0] = (char*)calloc(bufsize, sizeof(char));

    /* add a single flux variable y(i,j) forall i,j */
    for (int i = 0; i < nnodes; i++) {
        for (int j = 0; j < nnodes; j++) {
            if (i == j) continue;

            /* cost is zero for new variables, they matter for new constraints
             * only */
            snprintf(cname[0], bufsize, "y(%d)(%d)", i + 1, j + 1);
            double obj = 0;

            /* delete the constraints just adding this upper bound */
            if (j == 0) ub = 0;

            /* y_1i is actually = (n-1) x_ij so it can be at max n-1 */
            if (i == 0)
                ub = nnodes - 1;
            else
                ub = nnodes - 2;

            /* inject variable and test it's position (upos) inside CPLEX */
            if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &integer, cname)) {
                print_error("wrong CPXnewcols on x var.s");
            }
            if (CPXgetnumcols(env, lp) - 1 != ypos(i, j, nnodes)) {
                print_error(" wrong position for x var.s");
            }
        }
    }

    /* update the number of columns */
    inst->ncols += nnodes * nnodes;

    /* linking constraints */
    double rhs = 0.0;
    char sense = 'L';

    /* linking constraint: y_ij <= (n-2) x_ij forall i,j in V \ {1}
     * transformed in (2-n) x_ij + y_ij <= 0 */
    for (int i = 1; i < nnodes; i++) {
        for (int h = 1; h < nnodes; h++) {
            if (i == h) continue;

            /* fetch last row position in current model (better: the new one)
             * in CPLEX row == constraint */
            int lastrow = CPXgetnumrows(env, lp);
            snprintf(cname[0], bufsize, "xy(%d)(%d)", i + 1, h + 1);

            /* add the new constraint and change coeff from zero to 1 and 1-m */
            if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
                print_error("wrong CPXnewrows [degree]");
            }
            if (CPXchgcoef(env, lp, lastrow, xxpos(i, h, nnodes), 2 - nnodes)) {
                print_error("wrong CPXchgcoef [degree]");
            }
            if (CPXchgcoef(env, lp, lastrow, ypos(i, h, nnodes), 1.0)) {
                print_error("wrong CPXchgcoef [degree]");
            }
        }
    }

    /* node 1 out flow :
     * y_1j = (n-1) x_ij forall i in V \ {1}*/
    rhs = 0;
    sense = 'E';
    for (int i = 1; i < nnodes; i++) {
        snprintf(cname[0], bufsize, "y1%d", i + 1);
        int lastrow = CPXgetnumrows(env, lp);
        /* add the new constraint (with coefficent zero, null lhs) */
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
            print_error("wrong CPXnewrows [degree]");
        }

        if (CPXchgcoef(env, lp, lastrow, ypos(0, i, nnodes), 1.0)) {
            print_error("wrong CPXchgcoef [degree]");
        }

        if (CPXchgcoef(env, lp, lastrow, xpos(0, i + 1, nnodes),
                       1.0 - nnodes)) {
            print_error("wrong CPXchgcoef [degree]");
        }
    }

    /* node 1 in flow :
     * y_i1 = 0 forall i in V \ {1}*/
    /*
     *    rhs = 0;
     *    sense = 'E';
     *    for (int i = 1; i < nnodes; i++) {
     *        sprintf(cname[0], "y%d1", i + 1);
     *
     *        int lastrow = CPXgetnumrows(env, lp);
     *
     *        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
     *            print_error("wrong CPXnewrows [degree]");
     *        }
     *
     *        if (CPXchgcoef(env, lp, lastrow, ypos(i, 0, nnodes), 1.0)) {
     *            print_error("wrong CPXchgcoef [degree]");
     *        }
     *    }
     */
    /* alternative approach: sum i != 1 y_i1 = 0 */
    /*
     *    rhs = 0;
     *    sense = 'E';
     *    int lastrow = CPXgetnumrows(env, lp);
     *    sprintf(cname[0], "y1j");
     *
     *    [> add the new constraint (with coefficent zero, null lhs) <]
     *    if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
     *        print_error("wrong CPXnewrows [degree]");
     *    }
     *    for (int j = 1; j < nnodes; j++) {
     *        if (CPXchgcoef(env, lp, lastrow, ypos(j, 0, nnodes), 1.0)) {
     *            print_error("wrong CPXchgcoef [degree]");
     *        }
     *    }
     */
    /* alternative approach (used): set proper ub in variable! */

    /* flow control:
     * sum i != j  y_ij - sum h != j y_jh = 1 forall j in V \ {1}*/
    rhs = 1;
    sense = 'E';
    for (int j = 1; j < nnodes; j++) {
        int lastrow = CPXgetnumrows(env, lp);
        snprintf(cname[0], bufsize, "y(%d)", j + 1);

        /* add the new constraint (with coefficent zero, null lhs) */
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
            print_error("wrong CPXnewrows [degree]");
        }

        /* and change the coefficent if node i is adiacent to h */
        for (int i = 0; i < nnodes; i++) {
            /* the graph is complete so we only skip self loops */
            if (i == j) continue;

            /* change the coefficent from 0 to 1.0 */
            if (CPXchgcoef(env, lp, lastrow, ypos(i, j, nnodes), 1.0)) {
                print_error("wrong CPXchgcoef [degree]");
            }
        }

        /* and change the coefficent if node i is adiacent to h */
        for (int h = 0; h < nnodes; h++) {
            /* the graph is complete so we only skip self loops */
            if (h == j) continue;

            /* change the coefficent from 0 to 1.0 */
            if (CPXchgcoef(env, lp, lastrow, ypos(j, h, nnodes), -1.0)) {
                print_error("wrong CPXchgcoef [degree]");
            }
        }
    }
}

void add_GGlit_lazy_sec(CPXENVptr env, CPXLPptr lp, instance inst) {
    int nnodes = inst->nnodes;
    /* new integer variables
     * upper bound equal to m-1 */
    char integer = 'I';
    double lb = 0.0;
    double ub = nnodes - 1;

    /* ColumnNAME: array of array used to inject variables in CPLEX */
    char** cname = (char**)malloc(sizeof(char*));
    int bufsize = 100;
    cname[0] = (char*)calloc(bufsize, sizeof(char));

    /* add a single flux variable y(i,j) forall i,j */
    for (int i = 0; i < nnodes; i++) {
        for (int j = 0; j < nnodes; j++) {
            if (i == j) continue;

            /* cost is zero for new variables, they matter for new constraints
             * only */
            snprintf(cname[0], bufsize, "y(%d)(%d)", i + 1, j + 1);
            double obj = 0;

            /* inject variable and test it's position (upos) inside CPLEX */
            if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &integer, cname)) {
                print_error("wrong CPXnewcols on x var.s");
            }
            if (CPXgetnumcols(env, lp) - 1 != ypos(i, j, nnodes)) {
                print_error(" wrong position for x var.s");
            }
        }
    }

    /* update the number of columns */
    inst->ncols += nnodes * nnodes;

    int izero = 0;
    int index[2 * nnodes - 2];
    double value[2 * nnodes - 2];
    double rhs = 0.0;
    char sense = 'L';
    int nnz = 2;

    /* linking constraint: y_ij <= (n-1) x_ij forall i,j in V
     * transformed in (1-n) x_ij + y_ij <= 0 */
    for (int i = 0; i < nnodes; i++) {
        for (int j = 0; j < nnodes; j++) {
            if (i == j) continue;

            /* tighten constraint */
            int x_coeff;
            if (i != 0)
                x_coeff = 2 - nnodes;
            else
                x_coeff = 1 - nnodes;

            snprintf(cname[0], bufsize, "xy(%d)(%d)", i + 1, j + 1);

            /* build constraint equation */
            index[0] = xxpos(i, j, nnodes);
            value[0] = x_coeff;
            index[1] = ypos(i, j, nnodes);
            value[1] = 1.0;

            if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero,
                                      index, value, cname)) {
                print_error("wrong CPXlazyconstraints() for linking-GGlit");
            }
        }
    }

    /* sum j in V \ {1} y_1j = n - 1 */
    rhs = nnodes - 1;
    sense = 'E';
    nnz = nnodes - 1;
    izero = 0; /* index and value are setted starting from 1 */

    snprintf(cname[0], bufsize, "y1j");

    for (int j = 1; j < nnodes; j++) {
        index[j - 1] = ypos(0, j, nnodes);
        value[j - 1] = 1.0; /* this can be faster memsetted */
    }
    if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index,
                              value, cname)) {
        print_error("wrong CPXlazyconstraints() for linking-GGlit");
    }

    /* flow control:
     * sum i != j  y_ij - sum k != j y_jk = 1 forall j in V \ {1} */
    rhs = 1.0;
    sense = 'E';
    nnz = 2 * nnodes - 2;
    izero = 0;
    for (int j = 1; j < nnodes; j++) {
        snprintf(cname[0], bufsize, "y(%d)", j + 1);
        int u = 0;

        /* and change the coefficent if node i is adiacent to j */
        for (int i = 0; i < nnodes; i++) {
            if (i == j) continue;

            index[u] = ypos(i, j, nnodes);
            value[u] = 1.0;
            u++;
        }

        /* and change the coefficent if node j is adiacent to h */
        for (int k = 0; k < nnodes; k++) {
            if (k == j) continue;

            index[u] = ypos(j, k, nnodes);
            value[u] = -1.0;
            u++;
        }

        if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index,
                                  value, cname)) {
            print_error("wrong CPXlazyconstraints() for flow control");
        }
    }

    free(cname[0]);
    free(cname);
}

void perform_BENDERS(CPXENVptr env, CPXLPptr lp, instance inst, double* xstar,
                     struct timespec s, struct timespec e) {
    int nedges = inst->nnodes;

    /* create an adjacency list to track the edges .. */
    adjlist l = adjlist_create(nedges);
    /* .. and first fill it with the first solution found (with possible
     * subtours) */
    for (int i = 0; i < nedges; i++) {
        for (int j = i + 1; j < nedges; j++) {
            if (xstar[xpos(i, j, nedges)] > 0.5) {
                adjlist_add_edge(l, i, j);
            }
        }
    }

    /* iterate and add SEC until we obtain a singe tour */
    while (!adjlist_single_tour(l)) {
        /* add the constraints */
        add_BENDERS_sec(env, lp, l);

        /* some time has passed, update the timelimit the timelimit
         * we have to pass time in second! */
        long restime = inst->params->timelimit - stopwatch(&s, &e) / 1000.0;
        CPXsetdblparam(env, CPX_PARAM_TILIM, restime);

        /* solve again! */
        if (CPXmipopt(env, lp)) print_error("CPXmipopt() error");

        /* store the optimal solution found by CPLEX */
        if (CPXgetx(env, lp, xstar, 0, CPXgetnumcols(env, lp) - 1)) {
            print_error("CPXgetx() error\n");
        }

        /* completely reset the adjacency list without recreating the object and
         * fill it with the found solution */
        adjlist_hard_reset(l);
        for (int i = 0; i < nedges; i++) {
            for (int j = i + 1; j < nedges; j++) {
                if (xstar[xpos(i, j, nedges)] > 0.5) {
                    adjlist_add_edge(l, i, j);
                }
            }
        }
    }

    /* save the now complete model */
    char* filename;
    int bufsize = 100;
    filename = (char*)calloc(bufsize, sizeof(char));
    char* folder = model_folder_tostring(inst->model_folder);

    snprintf(filename, bufsize, "../data_%s/%s/%s.benders.lp", folder,
             inst->model_name, inst->model_name);

    CPXwriteprob(env, lp, filename, NULL);
    free(folder);
    free(filename);
}

void add_BENDERS_sec(CPXENVptr env, CPXLPptr lp, adjlist l) {
    /* add
     * sum i in V, j in V x_ij <= |V| - 1 for each V subtour of actual solution
     * to the model
     */
    int nnodes = l->N;

    char** cname = (char**)calloc(1, sizeof(char*));
    cname[0] = (char*)calloc(100, sizeof(char));

    /* iterate over all possible subtour nodes present in the solution and, for
     * a single subtour, add a SEC so no subtour will form during next
     * iterations */
    int* subtour;
    int subsize;
    int nsubtours = 0;
    while ((subtour = adjlist_get_subtour(l, &subsize)) != NULL) {
        nsubtours++;

        /* rhs of constraint is nodes in subtour -1, because of the other "one
         * edge entering and one exiting" for each node constrain, in a subset
         * of N nodes will be present N edges. Imposing N-1 fix this issue, at
         * least one edge has to exit the subset */
        const char sense = 'L';
        double rhs = (double)subsize - 1;

        /* fetch new row .. */
        int lastrow = CPXgetnumrows(env, lp);
        snprintf(cname[0], strlen(cname[0]), "benders_sec(%d)", subtour[0] + 1);
        /* .. and fix the sense and rhs */
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
            print_error("wrong CPXnewrows [benders]");
        }

        /* iterate over all possible ordered pair (symmetric model) of the
         * subtour and prevent that edge from forming*/
        for (int i = 0; i < subsize; i++) {
            for (int j = i + 1; j < subsize; j++) {
                if (CPXchgcoef(env, lp, lastrow,
                               xpos(subtour[i], subtour[j], nnodes), 1.0)) {
                    print_error("wrong CPXchgcoef [benders]");
                }
            }
        }

        free(subtour);
    }

    if (VERBOSE) printf("[VERBOSE] num subtour BENDERS %d\n", nsubtours);

    free(cname[0]);
    free(cname);
}

int CPXPUBLIC add_BENDERS_sec_callback_driver(CPXCALLBACKCONTEXTptr context,
                                              CPXLONG contextid,
                                              void* userhandle) {
    solution sol = (solution)userhandle;

    if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE)
        return add_BENDERS_sec_callback_candidate(context, sol);
    if (contextid == CPX_CALLBACKCONTEXT_RELAXATION)
        return add_BENDERS_sec_callback_relaxation(context, sol);

    print_error("cannot handle different contextids");
    return 1;
}

int CPXPUBLIC add_BENDERS_sec_callback_candidate(CPXCALLBACKCONTEXTptr context,
                                                 solution sol) {
    /* number of columns is n chooses 2 */
    int nedges = sol->nedges;
    int ncols = nedges * (nedges - 1) / 2;

    /* retreive actual solution */
    double* xstar = (double*)calloc(ncols, sizeof(double));
    double objval = CPX_INFBOUND;
    if (CPXcallbackgetcandidatepoint(context, xstar, 0, ncols - 1, &objval)) {
        print_error("CPXcallbackgetcandidatepoint error");
    }

    /* get node informations */
    int mynode = -1;
    int mythread = -1;
    double zbest;
    double incumbent = CPX_INFBOUND;
    CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODECOUNT, &mynode);
    CPXcallbackgetinfoint(context, CPXCALLBACKINFO_THREADID, &mythread);
    CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_SOL, &zbest);
    CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_SOL, &incumbent);
    if (VERBOSE && !SUPPRESS_CALLBACK) {
        printf("node information:\n");
        printf("- node index: %d\n", mynode);
        printf("- thread index: %d\n", mythread);
        printf("- zbest: %lf\n", zbest);
        printf("- incumbent: %lf\n", incumbent);
    }

    /* create an adjacency list to track the edges .. */
    adjlist l = adjlist_create(nedges);
    /* .. and first fill it with the first solution found (with possible
     * subtours) */
    for (int i = 0; i < nedges; i++) {
        for (int j = i + 1; j < nedges; j++) {
            if (xstar[xpos(i, j, nedges)] > 0.5) {
                adjlist_add_edge(l, i, j);
            }
        }
    }

    /* iterate over all possible subtour nodes present in the solution and, for
     * a single subtour, add a SEC so no subtour will form during next
     * iterations */
    int* subtour;
    int subsize;
    int nsubtours = 0;
    while ((subtour = adjlist_get_subtour(l, &subsize)) != NULL) {
        nsubtours++;

        /* rhs of constraint is nodes in subtour -1, because of the other "one
         * edge entering and one exiting" for each node constrain, in a subset
         * of N nodes will be present N edges. Imposing N-1 fix this issue, at
         * least one edge has to exit the subset */
        const char sense = 'L';
        double rhs = (double)subsize - 1;

        /* better to store indexes and positions of constraint's edges in
         * additional structures */
        int nnz = 0;
        const int izero = 0;
        int* index = (int*)malloc(ncols * sizeof(int));
        double* value = (double*)malloc(ncols * sizeof(double));

        /* iterate over all possible ordered pair (symmetric model) of the
         * subtour and prevent that edge from forming*/
        for (int i = 0; i < subsize; i++) {
            for (int j = i + 1; j < subsize; j++) {
                index[nnz] = xpos(subtour[i], subtour[j], nedges);
                value[nnz++] = 1.0;
            }
        }

        /* finally set the callback for rejecte the incumbent */
        if (rhs != nedges - 1) {
            if (CPXcallbackrejectcandidate(context, 1, nnz, &rhs, &sense,
                                           &izero, index, value)) {
                print_error("CPXcallbackrejectcandidate() error");
            }
        }

        free(index);
        free(value);
        free(subtour);
    }

    if (VERBOSE && !SUPPRESS_CALLBACK) {
        printf("[VERBOSE] num subtour BENDERS (callback) %d\n", nsubtours);
    }

    return 0;
}

int CPXPUBLIC add_BENDERS_sec_callback_relaxation(CPXCALLBACKCONTEXTptr context,
                                                  solution sol) {
    int nedges = sol->nedges;

    /* struct to pass info to doit_fn_callback */
    doit_fn_input data =
        (doit_fn_input)calloc(1, sizeof(struct doit_fn_input_t));
    data->nedges = nedges;
    data->context = context;

    /* number of columns is n chooses 2 */
    int ncols = nedges * (nedges - 1) / 2;

    /* get node informations */
    double* xstar = (double*)calloc(ncols, sizeof(double));
    double objval = CPX_INFBOUND;
    double epsilon = 0.1; /* set epsilon for CCcut_violated_cuts */
    if (CPXcallbackgetrelaxationpoint(context, xstar, 0, ncols - 1, &objval)) {
        print_error("CPXcallbackgetcandidatepoint error");
    }

    /* get node informations */
    int mynode = -1;
    int mythread = -1;
    double zbest;
    double incumbent = CPX_INFBOUND;
    CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODECOUNT, &mynode);
    CPXcallbackgetinfoint(context, CPXCALLBACKINFO_THREADID, &mythread);
    CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_SOL, &zbest);
    CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_SOL, &incumbent);
    if (VERBOSE && !SUPPRESS_CALLBACK) {
        printf("node information:\n");
        printf("- node index: %d\n", mynode);
        printf("- thread index: %d\n", mythread);
        printf("- zbest: %lf\n", zbest);
        printf("- incumbent: %lf\n", incumbent);
    }

    /* random discard the cut with 90% probability */
    // TODO(lugot): TUNE 90% hyperparameter
    if (mythread % 10 == 9) return 0;

    /* elist specify the vertices */
    int* elist = (int*)malloc(2 * ncols * sizeof(int));
    int loader = 0;
    for (int i = 0; i < nedges; i++) {
        for (int j = i + 1; j < nedges; j++) {
            elist[loader++] = i;
            elist[loader++] = j;
        }
    }

    /* componets infos */
    int ncomps = 0;
    int* comps = (int*)malloc(nedges * sizeof(int));
    int* compscount = (int*)malloc(nedges * sizeof(int));

    if (CCcut_connect_components(nedges, ncols, elist, xstar, &ncomps,
                                 &compscount, &comps)) {
        print_error("CCcut_connect_components error");
    }

    /* it seems it never happends */
    if (ncomps == 1) {
        if (VERBOSE) printf("[VERBOSE] add relaxation cut\n");
        if (CCcut_violated_cuts(nedges, ncols, elist, xstar, 2 - epsilon,
                                doit_fn_concorde, (void*)data))
            print_error("CCcut_violated_cuts error");
    }

    free(data);
    free(elist);
    free(xstar);
    free(comps);
    free(compscount);

    return 0;
}

int doit_fn_concorde(double cutval, int cutcount, int* cut, void* in) {
    doit_fn_input data = (doit_fn_input)in;
    double rhs = cutcount - 1.0;
    int nnz = 0;
    char sense = 'L';
    int purgeable = CPX_USECUT_FILTER;
    int local = 0;
    int izero = 0;

    double* value =
        (double*)calloc(cutcount * (cutcount - 1) / 2, sizeof(double));
    int* index = (int*)calloc(cutcount * (cutcount - 1) / 2, sizeof(int));

    // FIX: (data->sol)->nedges is not equal to the number of nodes!
    for (int i = 0; i < cutcount; i++) {
        for (int j = i + 1; j < cutcount; j++) {
            index[nnz] = xpos(cut[i], cut[j], data->nedges);
            value[nnz++] = 1.0;
        }
    }
    // TODO(magu): FIX segmentation fault: why?
    if (CPXcallbackaddusercuts(data->context, 1, nnz, &rhs, &sense, &izero,
                               index, value, &purgeable, &local)) {
        print_error("CPXcallbackaddusercuts() error");
    }
    // TODO(magu): TEST print cut

    free(index);
    free(value);
    return 0;
}

void perform_HARD_FIXING(CPXENVptr env, CPXLPptr lp, instance inst,
                         double* xstar, int perc) {
    assert(perc >= 0 && perc < 100);

    int nedges = inst->nnodes;
    int ncols = inst->ncols;

    /* TODO(lugot): FIX */
    srand(0);
    unsigned int seedp = 0L;

    /* objective tracking */
    double best_obj, prev_obj;
    CPXgetobjval(env, lp, &prev_obj);

    /* preparing some constants to quickly fix bounds */
    const char lbc = 'L';
    const char ubc = 'U';
    const double one = 1.0;
    const double zero = 0.0;

    /* perform HARD FIXING:
     * for each iteration fix randomly some edges, solve the model and release
     * the nodes fixed in order to warm start the next iteration */
    for (int iter = 0; iter < HF_ITERATIONS; iter++) {
        /* create an adjacency to track the fixed edges */
        adjlist l = adjlist_create(nedges);

        /* iterate over edges to randomly fix them
         * not a single for cause of arc tracking */
        for (int i = 0; i < nedges; i++) {
            for (int j = i + 1; j < nedges; j++) {
                int pos = xpos(i, j, nedges);

                /* check if the edge is an actual edge and select it */
                if (xstar[pos] > 0.5 && rand_r(&seedp) % 100 < perc) {
                    /* if (VERBOSE) { */
                    /*     printf("[VERBOSE] fixing arc (%d,%d)\n", i + 1, j +
                     * 1); */
                    /* } */

                    /* fix the edge by moving the lower bound to one */
                    CPXchgbds(env, lp, 1, &pos, &lbc, &one);

                    /* track which edges has been selected */
                    adjlist_add_edge(l, i, j);
                }
            }
        }

        /* iterate over tracked edges to retreive the loose ends:
         * let's create a simple SEC by fix that edge to zero */
        int i, j;
        while (adjlist_get_loose_ends(l, &i, &j)) {
            int pos = xpos(i, j, nedges);

            /* check if the edge is not part of a path (alone): in that case we
             * have to keep the lower bound to one to maintain consistency */
            double lb;
            CPXgetlb(env, lp, &lb, pos, pos);
            /* if lb is 1.0 means that we previously fix it -> skip */
            if (lb == 1.0) continue;

            /* if (VERBOSE) printf("[VERBOSE] removing (%d,%d)\n", i + 1, j +
             * 1); */
            CPXchgbds(env, lp, 1, &pos, &ubc, &zero);
        }

        /* save the model for analysis */
        if (EXTRA) CPXwriteprob(env, lp, "./hard_fixing", "LP");

        /* reset the timelimit: fraction of number of iterations */
        CPXsetdblparam(env, CPX_PARAM_TILIM,
                       inst->params->timelimit / HF_ITERATIONS);
        /* alternative: set the nodelimit */
        /* CPXsetintparam(env, CPX_PARAM_NODELIM, 100); */

        /* solve! and get solution */
        if (CPXmipopt(env, lp)) print_error("CPXmipopt() error");
        if (CPXgetx(env, lp, xstar, 0, ncols - 1)) {
            print_error("CPXgetx() error\n");
        }

        if (VERBOSE) {
            /* store the objective */
            CPXgetobjval(env, lp, &best_obj);

            printf(
                "[VERBOSE]: hard fixing status (iter %d out of %d)\n"
                "\tprev_obj: %.20lf\n"
                "\t act_obj: %.20lf\n",
                iter, HF_ITERATIONS, prev_obj, best_obj);

            prev_obj = best_obj;
        }

        if (EXTRA) {
            /* create a solution for the plot */
            solution sol = create_solution(inst, HARD_FIXING, nedges);
            sol->inst = inst;
            get_symmsol(xstar, sol);

            /* track which edges we fixed */
            int* edgecolors = (int*)calloc(nedges, sizeof(int));
            int u, v;
            adjlist_reset(l);
            while (adjlist_get_edge(l, &u, &v)) {
                for (int z = 0; z < nedges; z++) {
                    if (sol->edges[z].i == z && sol->edges[z].j == v) {
                        edgecolors[z] = 1;
                        break;
                    }
                }
            }

            /* save the plot */
            plot_graphviz(sol, edgecolors, iter);
        }

        /* free the adjlist, do it not for additional plot option */
        adjlist_free(l);

        /* relax the fixing: no need to check because most of the nodes are
         * fixed. Just checking if this is not the last iteration to provide a
         * coeherent model */
        if (iter != HF_ITERATIONS - 1) {
            for (int col = 0; col < ncols; col++) {
                CPXchgbds(env, lp, 1, &col, &lbc, &zero);
                CPXchgbds(env, lp, 1, &col, &ubc, &one);
            }
        }
    }

    if (VERBOSE) {
        printf(
            "[VERBOSE]: hard fixing status (on exit)\n"
            "\tprev_obj: %.20lf\n"
            "\t act_obj: %.20lf\n",
            prev_obj, best_obj);
    }
}

void perform_SOFT_FIXING(CPXENVptr env, CPXLPptr lp, instance inst,
                         double* xstar) {
    int nedges = inst->nnodes;
    int ncols = inst->ncols;

    double k = SF_INITIAL_K;
    double best_obj, prev_obj;
    CPXgetobjval(env, lp, &prev_obj);

    while (inst->params->timelimit > 0) {
        /* track the time for the current iteration */
        struct timespec s, e;
        s.tv_sec = e.tv_sec = -1;
        stopwatch(&s, &e);

        /* set the new timelimit and perform another resolution */
        CPXsetdblparam(env, CPX_PARAM_TILIM, inst->params->timelimit);

        char** cname = (char**)calloc(1, sizeof(char*));
        cname[0] = (char*)calloc(100, sizeof(char));

        // TODO(any): comment
        const char sense = 'L';
        const double rhs = k - nedges;

        /* fetch new row .. */
        int lastrow = CPXgetnumrows(env, lp);
        snprintf(cname[0], strlen(cname[0]), "soft_fixing");
        /* .. and fix the sense and rhs */
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
            print_error("wrong CPXnewrows [soft fixing]");
        }

        /* TODO(lugot): MODIFY double for for tracking fixed nodes */
        for (int i = 0; i < ncols; i++) {
            if (xstar[i] > 0.5) {
                /* change the coefficent from 0 to -1.0 for variables = 1 in
                 * xstar */
                if (CPXchgcoef(env, lp, lastrow, i, -1.0)) {
                    print_error("wrong CPXchgcoef [soft fixing]");
                }
            }
        }

        /* solve! and get solution */
        if (CPXmipopt(env, lp)) print_error("CPXmipopt() error");
        if (CPXgetx(env, lp, xstar, 0, ncols - 1)) {
            print_error("CPXgetx() error\n");
        }

        /* compute objective and compare it with the previous one to determine
         * if we need to enlarge the neighborhood size */
        CPXgetobjval(env, lp, &best_obj);

        if (VERBOSE) {
            printf(
                "[VERBOSE]: soft fixing status (bf k updating)\n"
                "\tremaining_time: %lf\n"
                "\tk: %lf\n"
                "\tprev_obj: %.20lf\n"
                "\t act_obj: %.20lf\n",
                inst->params->timelimit, k, prev_obj, best_obj);
        }

        /* if no improvment enlarge neigborhood size */
        if (fabs(prev_obj - best_obj) < EPSILON) k += SF_K_STEP;
        /* cut the computation if k too high */
        if (k >= SF_MAX_K) break;
        /* and update the objective */
        prev_obj = best_obj;

        /* update timelimit */
        inst->params->timelimit -= stopwatch(&s, &e) / 1000.0;

        if (inst->params->timelimit > 0)
            if (CPXdelrows(env, lp, lastrow, lastrow))
                print_error("CPXdelrows() error\n");
    }

    if (VERBOSE) {
        printf(
            "[VERBOSE]: soft fixing status (on exiting)\n"
            "\tremaining_time: %lf\n"
            "\tk: %lf\n"
            "\t act_obj: %.20lf\n",
            inst->params->timelimit, k, best_obj);
    }
}
