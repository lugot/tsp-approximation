#include "../include/solvers.h"

#include <assert.h>
#include <cplex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include "../include/adjlist.h"
#include "../include/globals.h"
#include "../include/model_builder.h"
#include "../include/tsp.h"
#include "../include/union_find.h"
#include "../include/utils.h"

void assert_correctness(solution sol);
void perform_BENDERS(CPXENVptr env, CPXLPptr lp, instance inst, solution sol,
                     int nedges, struct timeval start, struct timeval end,
                     double* xstar);
void perform_HARD_FIXING(CPXENVptr env, CPXLPptr lp, instance inst, int nedges,
                         int ncols, double* xstar, int perc);
void perform_SOFT_FIXING(CPXENVptr env, CPXLPptr lp, instance inst, int nedges,
                         int ncols, double* xstar);

solution TSPopt(instance inst, enum model_types model_type) {
    assert(inst != NULL);
    assert(inst->params != NULL && "no CPLEX params found");
    assert(inst->instance_type == TSP && "need TSP instance");

    /* create and populate solution */
    solution sol = create_solution(inst, model_type, inst->nnodes);
    int nedges = sol->nedges;
    sol->distance_time = compute_dist(inst);
    /* sol->timetype = inst->timetype; */

    /* open CPLEX model */
    int error;
    CPXENVptr env = CPXopenCPLEX(&error);
    CPXLPptr lp = CPXcreateprob(env, &error, "TSP");

    CPXgetdettime(env, &sol->start);

    /* set params:
     * - timelimit (timetype = 0 -> seconds, timetype = 1 -> ticks)
     * - EpInteger: CPLEX tollerance to declare a variable integer, important in
     * bigM
     * - EpRightHandSide: the less or equal satisfied up to this tollerance,
     * usually 1e-5, with bigM 1e-9 */

    // if(inst->timetype == 0)
    CPXsetdblparam(env, CPX_PARAM_TILIM, inst->params->timelimit);
    /*else
        CPXsetdblparam(env, CPX_PARAM_DETTILIM, inst->params->timelimit);*/

    // CPXsetintparam(env, CPX_PARAM_RANDOMSEED, seed); todo: add a seed to
    // Cplexrun (Passed by TSP opt?)
    CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0);
    CPXsetdblparam(env, CPX_PARAM_EPRHS, 1e-9);

    char* logfile;
    logfile = (char*)calloc(100, sizeof(char));
    if (inst->model_folder == TSPLIB)
        sprintf(logfile, "../data_tsplib/%s/%slog.txt", inst->model_name,
                inst->model_name);
    if (inst->model_folder == GENERATED)
        sprintf(logfile, "../data_generated/%s/%slog.txt", inst->model_name,
                inst->model_name);

    CPXsetlogfilename(env, logfile, "w");

    /* populate enviorment with model data */
    sol->build_time = build_tsp_model(env, lp, inst, model_type);

    /* model preprocessing: add some callbacks or set params before execution */
    switch (model_type) {
        case BENDERS_CALLBACK: {
            CPXLONG contextid =
                CPX_CALLBACKCONTEXT_CANDIDATE | CPX_CALLBACKCONTEXT_RELAXATION;
            if (CPXcallbacksetfunc(env, lp, contextid,
                                   add_BENDERS_sec_callback_driver, sol)) {
                print_error("CPXcallbacksetfunc() error");
            }
            break;
        }

        case HARD_FIXING: {
            CPXsetintparam(env, CPX_PARAM_NODELIM, 10);

            // if(inst->timetype == 0)
            CPXsetdblparam(env, CPX_PARAM_TILIM, inst->params->timelimit / 10);
            /*else
                    CPXsetdblparam(env, CPX_PARAM_DETTILIM,
               inst->params->timelimit/10);*/

            CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE;
            if (CPXcallbacksetfunc(env, lp, contextid,
                                   add_BENDERS_sec_callback_driver, sol)) {
                print_error("CPXcallbacksetfunc() error");
            }
            break;
        }

        default:
            /* no preprocessing steps */
            break;
    }

    struct timeval start, end;
    gettimeofday(&start, NULL);
    gettimeofday(&end, NULL); /* warning supperssor: avoid uninitialized var */

    /* solve! */
    if (CPXmipopt(env, lp)) print_error("CPXmipopt() error");

    /* store the optimal solution found by CPLEX */
    int ncols = CPXgetnumcols(env, lp);
    double* xstar = (double*)calloc(ncols, sizeof(double));
    if (CPXgetx(env, lp, xstar, 0, ncols - 1)) print_error("CPXgetx() error\n");

    /* retreive solution or perform future steps */
    switch (model_type) {
        case NOSEC:
        case BENDERS_CALLBACK:
            get_symmsol(xstar, sol);
            break;

        case BENDERS:
            perform_BENDERS(env, lp, inst, sol, nedges, start, end, xstar);
            get_symmsol(xstar, sol);
            break;

        case HARD_FIXING:
            perform_HARD_FIXING(env, lp, inst, nedges, ncols, xstar, 80);
            get_symmsol(xstar, sol);
            break;

        case SOFT_FIXING:
            perform_SOFT_FIXING(env, lp, inst, nedges, ncols, xstar);
            get_symmsol(xstar, sol);
            break;

        case MTZ_STATIC:
        case MTZ_LAZY:
        case GGLIT_STATIC:
        case GGLECT_STATIC:
        case GGLIT_LAZY:
            get_asymmsol(xstar, sol);
            break;

        case OPTIMAL_TOUR:
            assert(model_type != OPTIMAL_TOUR &&
                   "tried to solve an optimal tour instance");
            break;
    }

    free(xstar);

    CPXgetdettime(env, &sol->end);

    gettimeofday(&end, NULL);

    // if(inst->timetype == 0)
    sol->solve_time = (end.tv_sec - start.tv_sec);  // sec
    /*else
            sol->solve_time = sol->end - sol->start; //ticks*/

    /* retreive the min cost */
    CPXgetobjval(env, lp, &sol->zstar);
    /*sol->zstar = compute_zstar(inst, sol);*/

    /* assert_correctness(sol); // TODO(any): waiting for rewriting */

    /* add the solution to the pool associated with it's instance */
    add_solution(inst, sol);

    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);

    return sol;
}

void perform_BENDERS(CPXENVptr env, CPXLPptr lp, instance inst, solution sol,
                     int nedges, struct timeval start, struct timeval end,
                     double* xstar) {
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

        // TODO(any): time management
        /* some time has passed, updat the timelimit */
        gettimeofday(&end, NULL);
        /*if(inst->timetype == 0)
        {*/
        int restime = inst->params->timelimit - (end.tv_sec - start.tv_sec);
        CPXsetdblparam(env, CPX_PARAM_TILIM, restime);
        /*}
        else
        {
                int restime =
                inst->params->timelimit - (end.tv_sec -
        start.tv_sec)*TICKS_PER_SECOND; CPXsetdblparam(env,
        CPX_PARAM_DETTILIM, restime);
        }*/

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
    char* filename; //TODO (any): rewrite
    filename = (char*)calloc(100, sizeof(char));
    if (inst->model_folder == TSPLIB)
        snprintf(filename, 28 + 2 * strlen(inst->model_name),
                 "../data_tsplib/%s/%s.benders.lp", inst->model_name,
                 inst->model_name);
    if (inst->model_folder == GENERATED)
        snprintf(filename, 31 + 2 * strlen(inst->model_name),
                 "../data_generated/%s/%s.benders.lp", inst->model_name,
                 inst->model_name);
    CPXwriteprob(env, lp, filename, NULL);
    free(filename);
}

void perform_HARD_FIXING(CPXENVptr env, CPXLPptr lp, instance inst, int nedges,
                         int ncols, double* xstar, int perc) {
    assert(perc >= 0 && perc < 100);

    /* TODO(lugot): FIX */
    srand(0);
    unsigned int seedp = 0L;

    /* preparing some constants to quickly fix bounds */
    const char lbc = 'L';
    const char ubc = 'U';
    const double one = 1.0;
    const double zero = 0.0;

    /* perform HARD FIXING:
     * for each iteration fix randomly some edges, solve the model and release
     * the nodes fixed in order to warm start the next iteration */
    int iterations = 10;
    for (int iter = 0; iter < iterations; iter++) {
        if (VERBOSE) {
            printf("[VERBOSE] hard fixing iter %d out of %d\n", iter,
                   iterations);
        }

        /* create an adjacency to track the fixed edges */
        adjlist l = adjlist_create(nedges);

        /* iterate over edges to randomly fix them
         * not a single for cause of arc tracking */
        for (int i = 0; i < nedges; i++) {
            for (int j = i + 1; j < nedges; j++) {
                int pos = xpos(i, j, nedges);

                /* check if the edge is an actual edge and select it */
                if (xstar[pos] > 0.5 && rand_r(&seedp) % 100 < perc) {
                    if (VERBOSE) {
                        printf("[VERBOSE] fixing arc (%d,%d)\n", i + 1, j + 1);
                    }

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

            if (VERBOSE) printf("[VERBOSE] removing (%d,%d)\n", i + 1, j + 1);
            CPXchgbds(env, lp, 1, &pos, &ubc, &zero);
        }

        if (EXTRA) CPXwriteprob(env, lp, "./hard_fixing", "LP");

        // TODO(any): fix time
        CPXsetdblparam(env, CPX_PARAM_TILIM,
                       inst->params->timelimit / iterations);
        /*else
                CPXsetdblparam(env, CPX_PARAM_DETTILIM,
           inst->params->timelimit/repeat);*/

        CPXsetintparam(env, CPX_PARAM_NODELIM, 100);

        /* solve! and get solution */
        if (CPXmipopt(env, lp)) print_error("CPXmipopt() error");
        if (CPXgetx(env, lp, xstar, 0, ncols - 1)) {
            print_error("CPXgetx() error\n");
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
                for (int i = 0; i < nedges; i++) {
                    if (sol->edges[i].i == u && sol->edges[i].j == v) {
                        edgecolors[i] = 1;
                        break;
                    }
                }
            }

            /* save the plot */
            /* plot_solution_graphviz(sol, edgecolors, iter); TODO */
        }

        /* free the adjlist, do it not for additional plot option */
        adjlist_free(l);

        /* relax the fixing: no need to check because most of the nodes are
         * fixed. Just checking if this is not the last iteration to provide a
         * coeherent model */
        if (iter != iterations - 1) {
            for (int col = 0; col < ncols; col++) {
                CPXchgbds(env, lp, 1, &col, &lbc, &zero);
                CPXchgbds(env, lp, 1, &col, &ubc, &one);
            }
        }
    }
}

void perform_SOFT_FIXING(CPXENVptr env, CPXLPptr lp, instance inst, int nedges,
                         int ncols, double* xstar) {
    double timeremaining = 9 * inst->params->timelimit / 10;

    CPXsetdblparam(env, CPX_PARAM_NODELIM, CPX_INFBOUND);

    double k = 5.0;
    int noImproveFlag = 0;

    double bestVal;

    while (timeremaining > 0) {
        struct timeval start, end;
        gettimeofday(&start, NULL);

        double pvsObjVal = bestVal;
        CPXgetobjval(env, lp, &bestVal);

        if (pvsObjVal == bestVal)
            noImproveFlag = 1;
        else
            noImproveFlag = 0;

        printf("Soft fixing %f \n", timeremaining);
        CPXsetdblparam(env, CPX_PARAM_TILIM, timeremaining);

        if (noImproveFlag == 1)
            k += 2;
        else
            k = 5.0;

        char sense = 'L';
        double rhs = k - nedges;

        /* ColumnNAME: array of array used to inject variables in CPLEX */
        char** cname = (char**)calloc(1, sizeof(char*));
        cname[0] = (char*)calloc(100, sizeof(char));

        int lastrow = CPXgetnumrows(env, lp);
        snprintf(cname[0], strlen(cname[0]), "soft_fixing");

        /* add the new constraint (with coefficent zero, null lhs) */
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
            print_error("wrong CPXnewrows [degree]");
        }

        for (int i = 0; i < ncols; i++) {
            if (xstar[i] < 0.5) {
                /* change the coefficent from 0 to 1.0 for variables = 0 in
                 * xstar */
                if (CPXchgcoef(env, lp, lastrow, i, 1.0)) {
                    print_error("wrong CPXchgcoef [degree]");
                }
            } else {
                /* change the coefficent from 0 to -1.0 for variables = 1 in
                 * xstar */
                if (CPXchgcoef(env, lp, lastrow, i, -1.0)) {
                    print_error("wrong CPXchgcoef [degree]");
                }
            }
        }

        if (CPXmipopt(env, lp)) print_error("CPXmipopt() error");

        if (CPXgetx(env, lp, xstar, 0, ncols - 1))
            print_error("CPXgetx() error\n");

        gettimeofday(&end, NULL);
        double seconds = (double)(end.tv_sec - start.tv_sec);
        timeremaining -= seconds;

        if (timeremaining > 0)
            if (CPXdelrows(env, lp, lastrow, lastrow))
                print_error("CPXdelrows() error\n");
    }
}

void get_symmsol(double* xstar, solution sol) {
    /* index over selected edges */
    int k = 0;

    /* check for edges j>i if is selected */
    for (int i = 0; i < sol->nedges; i++) {
        for (int j = i + 1; j < sol->nedges; j++) {
            if (xstar[xpos(i, j, sol->nedges)] > 0.5) {
                /* then store it in the solution*/
                sol->edges[k++] = (edge){i, j};
            }
        }
    }

    // TODO(any): check this
    /*assert(k == nedges && "not enought edges CPLEX solution");*/
}

void get_asymmsol(double* xstar, solution sol) {
    /* index over selected edges */
    int k = 0;

    /* check for edges i, j if is selected */
    for (int i = 0; i < sol->nedges; i++) {
        for (int j = 0; j < sol->nedges; j++) {
            if (xstar[xpos(i, j, sol->nedges)] > 0.5) {
                /* then store it in the solution*/
                sol->edges[k++] = (edge){i, j};
            }
        }
    }

    // TODO(any): check this
    /*assert(k == nedges && "not enought edges CPLEX solution");*/
}

/* void assert_correctness(solution sol) { */
/*     TODO(any): rewrite why not */
/*     int visited = 0; */
/*     int act = 0; */
/*     do { */
/*         visited++; */
/*     } while ((act = sol->link[act]) != 0); */

/*     assert(visited == sol->nedges); */
/* } */
