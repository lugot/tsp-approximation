#include "../include/models/gg.h"

#include <cplex.h>

#include "../include/utils.h"

void add_GG_variables(CPXENVptr env, CPXLPptr lp, instance inst) {
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
}

void add_GGlit_static_sec(CPXENVptr env, CPXLPptr lp, instance inst) {
    int nnodes = inst->nnodes;

    /* ColumnNAME: array of array used to inject variables in CPLEX */
    char** cname = (char**)malloc(sizeof(char*));
    int bufsize = 100;
    cname[0] = (char*)calloc(bufsize, sizeof(char));

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

    char** cname = (char**)calloc(1, sizeof(char*));
    int bufsize = 100;
    cname[0] = (char*)calloc(bufsize, sizeof(char));

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

    free(cname[0]);
    free(cname);
}

void add_GGlit_lazy_sec(CPXENVptr env, CPXLPptr lp, instance inst) {
    int nnodes = inst->nnodes;

    /* ColumnNAME: array of array used to inject variables in CPLEX */
    char** cname = (char**)malloc(sizeof(char*));
    int bufsize = 100;
    cname[0] = (char*)calloc(bufsize, sizeof(char));

    int izero = 0;
    int* index;
    double* value;
    double rhs = 0.0;
    char sense = 'L';
    int nnz = 2;

    index = (int*)malloc((2 * nnodes - 2) * sizeof(int));
    value = (double*)malloc((2 * nnodes - 2) * sizeof(double));

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

void add_GGlect_lazy_sec(CPXENVptr env, CPXLPptr lp, instance inst) {
    int nnodes = inst->nnodes;

    /* ColumnNAME: array of array used to inject variables in CPLEX */
    char** cname = (char**)calloc(1, sizeof(char*));
    int bufsize = 100;
    cname[0] = (char*)calloc(bufsize, sizeof(char));

    /* linking constraints */
    int izero = 0;
    double rhs = 0.0;
    char sense = 'L';
    int* index;
    double* value;
    int nnz = 2;

    index = (int*)malloc((2 * nnodes - 2) * sizeof(int));
    value = (double*)malloc((2 * nnodes - 2) * sizeof(double));

    /* linking constraint: y_ij <= (n-2) x_ij forall i,j in V \ {1}
     * transformed in (2-n) x_ij + y_ij <= 0 */
    for (int i = 1; i < nnodes; i++) {
        for (int h = 1; h < nnodes; h++) {
            if (i == h) continue;

            /* fetch last row position in current model (better: the new one)
             * in CPLEX row == constraint */
            snprintf(cname[0], bufsize, "xy(%d)(%d)", i + 1, h + 1);

            index[0] = ypos(i, h, nnodes);
            value[0] = 1;
            index[1] = xxpos(i, h, nnodes);
            value[1] = 2 - nnodes;

            if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero,
                                      index, value, cname)) {
                print_error("wrong CPXlazyconstraints() for linking-GGlit");
            }
        }
    }

    /* node 1 out flow :
     * y_1j = (n-1) x_ij forall i in V \ {1}*/
    rhs = 0;
    sense = 'E';
    for (int j = 1; j < nnodes; j++) {
        snprintf(cname[0], bufsize, "y1%d", j + 1);
        /* add the new constraint (with coefficent zero, null lhs) */

        index[0] = ypos(0, j, nnodes);
        value[0] = 1;
        index[1] = xxpos(0, j, nnodes);
        value[1] = 1 - nnodes;

        if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index,
                                  value, cname)) {
            print_error("wrong CPXlazyconstraints() for linking-GGlit");
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
