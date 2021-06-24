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
#include "../include/models/mtz.h"
#include "../include/models/gg.h"
#include "../include/solvers.h"
#include "../include/tsp.h"
#include "../include/union_find.h"
#include "../include/utils.h"

void add_symm_variables(CPXENVptr env, CPXLPptr lp, instance inst);
void add_asymm_variables(CPXENVptr env, CPXLPptr lp, instance inst);

void add_symm_constraints(CPXENVptr env, CPXLPptr lp, instance inst);
void add_asymm_constraints(CPXENVptr env, CPXLPptr lp, instance inst);

void add_deg2_sec(CPXENVptr env, CPXLPptr lp, instance inst, modes mode);
void add_deg3_sec(CPXENVptr env, CPXLPptr lp, instance inst, modes mode);

/* void add_GGlit_static_deg2_sec(CPXENVptr env, CPXLPptr lp, instance inst); */
/*  */
/* int CPXPUBLIC add_BENDERS_sec_callback_candidate(CPXCALLBACKCONTEXTptr
 * context, */
/*                                                  solution sol); */
/* int CPXPUBLIC add_BENDERS_sec_callback_relaxation(CPXCALLBACKCONTEXTptr
 * context, */
/*                                                   solution sol); */
/* int doit_fn_concorde(double cutval, int cutcount, int* cut, void* in); */

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
            add_MTZ_variables(env, lp, inst);
            add_asymm_constraints(env, lp, inst);
            add_MTZ_static_sec(env, lp, inst);
            add_deg2_sec(env, lp, inst, STATIC);
            break;
        case MTZ_LAZY:
            add_asymm_variables(env, lp, inst);
            add_MTZ_variables(env, lp, inst);
            add_asymm_constraints(env, lp, inst);
            add_MTZ_lazy_sec(env, lp, inst);
            break;
        case MTZ_LAZY_DEG2:
            add_asymm_variables(env, lp, inst);
            add_MTZ_variables(env, lp, inst);
            add_asymm_constraints(env, lp, inst);
            add_MTZ_lazy_sec(env, lp, inst);
            add_deg2_sec(env, lp, inst, LAZY);
            break;
        case MTZ_LAZY_DEG3:
            add_asymm_variables(env, lp, inst);
            add_MTZ_variables(env, lp, inst);
            add_asymm_constraints(env, lp, inst);
            add_MTZ_lazy_sec(env, lp, inst);
            add_deg2_sec(env, lp, inst, LAZY);
            add_deg3_sec(env, lp, inst, LAZY);
            break;
        case MTZ_INDICATOR:
            add_asymm_variables(env, lp, inst);
            add_MTZ_variables(env, lp, inst);
            add_asymm_constraints(env, lp, inst);
            add_MTZ_indicator_sec(env, lp, inst);
            break;
        case GGLIT_STATIC:
            add_asymm_variables(env, lp, inst);
            add_GG_variables(env, lp, inst);
            add_asymm_constraints(env, lp, inst);
            add_GGlit_static_sec(env, lp, inst);
            break;
        case GGLECT_STATIC:
            add_asymm_variables(env, lp, inst);
            add_GG_variables(env, lp, inst);
            add_asymm_constraints(env, lp, inst);
            add_GGlect_static_sec(env, lp, inst);
            break;
        case GGLIT_LAZY:
            add_asymm_variables(env, lp, inst);
            add_GG_variables(env, lp, inst);
            add_asymm_constraints(env, lp, inst);
            add_GGlit_lazy_sec(env, lp, inst);
            break;
        case GGLECT_LAZY:
            add_asymm_variables(env, lp, inst);
            add_GG_variables(env, lp, inst);
            add_asymm_constraints(env, lp, inst);
            add_GGlect_lazy_sec(env, lp, inst);
            break;
        case GGLIT_STATIC_DEG2:
            add_asymm_variables(env, lp, inst);
            add_GG_variables(env, lp, inst);
            add_asymm_constraints(env, lp, inst);
            add_GGlit_static_sec(env, lp, inst);
            add_deg2_sec(env, lp, inst, STATIC);
            break;

        case MST:
        case GRASP:
        case GREEDY:
        case EXTRA_MILEAGE:
        case OPTIMAL_TOUR:
        case VNS:
        case TABU_SEACH:
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

void add_deg2_sec(CPXENVptr env, CPXLPptr lp, instance inst, modes mode) {
    int nnodes = inst->nnodes;
    /* add lazy constraints 1.0 * x_ij + 1.0 * x_ji <= 1
     * for each arc (i,j) with i < j */

    char** cname = (char**)calloc(1, sizeof(char*));
    cname[0] = (char*)calloc(100, sizeof(char));

    int izero = 0;
    int index[2];
    double value[2];
    double rhs = 1.0;
    char sense = 'L';
    int nnz = 2;

    for (int i = 0; i < nnodes; i++) {
        for (int j = i + 1; j < nnodes; j++) {
            snprintf(cname[0], strlen(cname[0]), "SEC on node pair (%d,%d)",
                     i + 1, j + 1);

            /* build constraint equation */
            index[0] = xxpos(i, j, nnodes);
            value[0] = 1.0;
            index[1] = xxpos(j, i, nnodes);
            value[1] = 1.0;

            switch (mode) {
                case STATIC: {
                    if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero,
                                   index, value, NULL, cname)) {
                        print_error("wrong CPXaddrows() for deg2 SEC");
                    }
                    break;
                }
                case LAZY: {
                    if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense,
                                              &izero, index, value, cname)) {
                        print_error("wrong CPXlazyconstraints on 2-node SECs");
                    }
                    break;
                }
            }
        }
    }
}
void add_deg3_sec(CPXENVptr env, CPXLPptr lp, instance inst, modes mode) {
    int nnodes = inst->nnodes;
    /* add lazy constraints x_ij + x_jk + x_ki <= 2
     * for each i, j, k with i < j < k and also on the opposite path */

    char** cname = (char**)calloc(1, sizeof(char*));
    cname[0] = (char*)calloc(100, sizeof(char));

    int izero = 0;
    int index[3];
    double value[3];
    double rhs = 2.0;
    char sense = 'L';
    int nnz = 3;

    for (int i = 0; i < nnodes; i++) {
        for (int j = i + 1; j < nnodes; j++) {
            for (int k = j + 1; k < nnodes; k++) {
                snprintf(cname[0], strlen(cname[0]),
                         "SEC on node triplet (%d,%d,%d)", i + 1, j + 1, k + 1);

                /* build constraint equation */
                index[0] = xxpos(i, j, nnodes);
                value[0] = 1.0;
                index[1] = xxpos(j, k, nnodes);
                value[1] = 1.0;
                index[2] = xxpos(k, i, nnodes);
                value[2] = 1.0;
                switch (mode) {
                    case STATIC: {
                        if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero,
                                       index, value, NULL, cname)) {
                            print_error("wrong CPXaddrows() for deg2 SEC");
                        }
                        break;
                    }
                    case LAZY: {
                        if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense,
                                                  &izero, index, value,
                                                  cname)) {
                            print_error(
                                "wrong CPXlazyconstraints on 3-node SECs");
                        }
                        break;
                    }
                }

                /* build constraint equation, inverse */
                index[0] = xxpos(i, k, nnodes);
                value[0] = 1.0;
                index[1] = xxpos(k, j, nnodes);
                value[1] = 1.0;
                index[2] = xxpos(j, i, nnodes);
                value[2] = 1.0;
                switch (mode) {
                    case STATIC: {
                        if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero,
                                       index, value, NULL, cname)) {
                            print_error("wrong CPXaddrows() for deg2 SEC");
                        }
                        break;
                    }
                    case LAZY: {
                        if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense,
                                                  &izero, index, value,
                                                  cname)) {
                            print_error(
                                "wrong CPXlazyconstraints on 3-node SECs");
                        }
                        break;
                    }
                }
            }
        }
    }
}
