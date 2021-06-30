#include "../../include/models/mtz.h"

#include <string.h>

#include "../../include/utils.h"

void add_MTZ_variables(CPXENVptr env, CPXLPptr lp, instance inst) {
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
}

void add_MTZ_static_sec(CPXENVptr env, CPXLPptr lp, instance inst) {
    int nnodes = inst->nnodes;

    /* ColumnNAME: array of array used to inject variables in CPLEX */
    char** cname = (char**)calloc(1, sizeof(char*));
    cname[0] = (char*)calloc(100, sizeof(char));

    int izero = 0;
    int index[3];
    double value[3];
    double big_M = nnodes - 1.0;
    double rhs = big_M - 1.0;
    char sense = 'L';
    int nnz = 3;

    /* add static constraints  1.0 * u_i - 1.0 * u_j + M * x_ij <= M - 1
     * for each arc (i,j) not touching node 0 */
    for (int i = 1; i < nnodes; i++) {
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
                print_error("wrong CPXaddrows() for u-consistency");
            }
        }
    }

    free(cname[0]);
    free(cname);
}

void add_MTZ_lazy_sec(CPXENVptr env, CPXLPptr lp, instance inst) {
    int nnodes = inst->nnodes;

    /* ColumnNAME: array of array used to inject variables in CPLEX */
    char** cname = (char**)calloc(1, sizeof(char*));
    cname[0] = (char*)calloc(100, sizeof(char));

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

    free(cname[0]);
    free(cname);
}

void add_MTZ_indicator_sec(CPXENVptr env, CPXLPptr lp, instance inst) {
    int nnodes = inst->nnodes;

    /* ColumnNAME: array of array used to inject variables in CPLEX */
    char** cname = (char**)calloc(1, sizeof(char*));
    cname[0] = (char*)calloc(100, sizeof(char));

    /* update the number of columns */
    inst->ncols += nnodes;

    int index[2];
    double value[2];
    double rhs = -1.0;
    char sense = 'L';
    int nnz = 2;

    /* add static constraints  1.0 * u_i - 1.0 * u_j + M * x_ij <= M - 1
     * for each arc (i,j) not touching node 0 */
    for (int i = 1; i < nnodes; i++) {
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

            if (CPXaddindconstr(env, lp, xxpos(i, j, nnodes), 0, nnz, rhs,
                                sense, index, value, cname[0])) {
                print_error("wrong CPXaddindconstr() for u-consistency");
            }
        }
    }

    free(cname[0]);
    free(cname);
}
