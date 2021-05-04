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

    /* add callback if required */
    if (model_type == BENDERS_CALLBACK) {
        CPXLONG contextid =
            CPX_CALLBACKCONTEXT_CANDIDATE | CPX_CALLBACKCONTEXT_RELAXATION;
        if (CPXcallbacksetfunc(env, lp, contextid,
                               add_BENDERS_sec_callback_driver, sol)) {
            print_error("CPXcallbacksetfunc() error");
        }
    }

    if (model_type == HARD_FIXING) {
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

    switch (model_type) {
        case NOSEC:
        case BENDERS_CALLBACK:
            get_symmsol(xstar, nedges, sol->edges, sol->link);
            break;

        case BENDERS:
            get_symmsol(xstar, nedges, sol->edges, sol->link);
            perform_BENDERS(env, lp, inst, sol, nedges, start, end, xstar);
            /*get_symmsol(xstar, nedges, sol->edges, sol->link); useless */
            break;

        case HARD_FIXING:
            perform_HARD_FIXING(env, lp, inst, nedges, ncols, xstar);
            get_symmsol(xstar, nedges, sol->edges, sol->link);
            break;

        case MTZ_STATIC:
        case MTZ_LAZY:
        case GGLIT_STATIC:
        case GGLECT_STATIC:
        case GGLIT_LAZY:
            get_asymmsol(xstar, nedges, sol->edges, sol->link);
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

    /* assert single cycle */
    /*print_solution(sol, 1);*/

    assert_correctness(sol);

    /* add the solution to the pool associated with it's instance */
    add_solution(inst, sol);

    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);

    return sol;
}

void perform_BENDERS(CPXENVptr env, CPXLPptr lp, instance inst, solution sol,
                     int nedges, struct timeval start, struct timeval end,
                     double* xstar) {
    while (!visitable(sol->link, nedges)) {
        /* add the constraints */
        add_BENDERS_sec(env, lp, sol);

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

        /* solve! */
        if (CPXmipopt(env, lp)) print_error("CPXmipopt() error");

        /* store the optimal solution found by CPLEX */
        int num_cols = CPXgetnumcols(env, lp);

        memset(xstar, '\0', num_cols * sizeof(double));
        if (CPXgetx(env, lp, xstar, 0, num_cols - 1))
            print_error("CPXgetx() error\n");

        /* retrive the solution */
        get_symmsol(xstar, nedges, sol->edges, sol->link);
    }

    /* save the complete model */
    char* filename;
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
                         int ncols, double* xstar) {
    srand(0);
    unsigned int seedp;
    int repeat = 10;

    char lbc, ubc;
    lbc = 'L';
    ubc = 'U';
    int pos;
    double one = 1.0;
    double zero = 0.0;

    int bsize = 100;
    char* version = (char*)calloc(bsize, sizeof(char));

    for (int iter = 0; iter < repeat; iter++) {
        adjlist l = adjlist_create(nedges);

        /* fix nodes */
        for (int i = 0; i < nedges; i++) {
            for (int j = i + 1; j < nedges; j++) {
                if (rand_r(&seedp) % 10 > 5 &&
                    xstar[xpos(i, j, nedges)] > 1.0 - 1e-9) {
                    if (VERBOSE)
                        printf("[Verbose] fixing (%d,%d)\n", i + 1, j + 1);

                    pos = xpos(i, j, nedges);
                    CPXchgbds(env, lp, 1, &pos, &lbc, &one);

                    adjlist_add_arc(l, i, j);
                }
            }
        }

        int i, j;
        while (adjlist_loose_ends(l, &i, &j)) {
            if (VERBOSE) printf("[Verbose] removing (%d,%d)\n", i + 1, j + 1);

            pos = xpos(i, j, nedges);
            /* CPXchgbds(env, lp, 1, &pos, &ubc, &zero); */
        }

        CPXwriteprob(env, lp, "./wow", "LP");
        // if(inst->timetype == 0)
        CPXsetdblparam(env, CPX_PARAM_TILIM, inst->params->timelimit / repeat);
        /*else
                CPXsetdblparam(env, CPX_PARAM_DETTILIM,
           inst->params->timelimit/repeat);*/

        CPXsetintparam(env, CPX_PARAM_NODELIM, 100);

        if (CPXmipopt(env, lp)) print_error("CPXmipopt() error");

        if (CPXgetx(env, lp, xstar, 0, ncols - 1)) {
            print_error("CPXgetx() error\n");
        }

        /* create a solution for the plot */
        solution sol = create_solution(inst, HARD_FIXING, nedges);
        sol->inst = inst;
        get_symmsol(xstar, nedges, sol->edges, sol->link);

        /* track which edges we fixed */
        int* edgecolors = (int*)calloc(nedges, sizeof(int));
        int u, v;
        adjlist_reset(l);
        while (adjlist_get_arc(l, &u, &v)) {
            for (int i = 0; i < nedges; i++) {
                if (sol->edges[i].i == u && sol->edges[i].j == v) {
                    edgecolors[i] = 1;
                    break;
                }
            }
        }

        /* save the plot */
        printf("diocane%d\n", iter);
        if (iter == 1) plot_solution_graphviz(sol, edgecolors, iter);

        adjlist_free(l);

        /* free nodes: no check because most of the nodes are fixed */
        if (iter != repeat - 1) {
            for (int j = 0; j < ncols; j++) {
                CPXchgbds(env, lp, 1, &j, &lbc, &zero);
                CPXchgbds(env, lp, 1, &j, &ubc, &one);
            }
        }
    }

    free(version);
}

void get_symmsol(double* xstar, int nedges, edge* edges, int* link) {
    assert(link != NULL);

    /* set up parent as linked list, no check on size */
    for (int i = 0; i < nedges; i++) link[i] = i;

    /* index over selected edges */
    int k = 0;

    /* check for edges j>i if is selected */
    for (int i = 0; i < nedges; i++) {
        for (int j = i + 1; j < nedges; j++) {
            if (xstar[xpos(i, j, nedges)] > 0.5) {
                /* link is mandatory, edges is optional */
                if (edges != NULL) edges[k] = (edge){i, j};
                k++;

                /* actually a shifted linked list */
                if (!reachable(link, i, j) && !reachable(link, j, i))
                    swap(&link[i], &link[j]);
            }
        }
    }

    /*assert(k == nedges && "not enought edges CPLEX solution");*/
}

void get_asymmsol(double* xstar, int nedges, edge* edges, int* link) {
    assert(link != NULL);

    /* set up parent as linked list */
    for (int i = 0; i < nedges; i++) link[i] = i;

    /* index over selected edges */
    int k = 0;

    /* check for all edges if is selected */
    for (int i = 0; i < nedges; i++) {
        for (int j = 0; j < nedges; j++) {
            if (xstar[xxpos(i, j, nedges)] > 0.5) {
                /* link is mandatory, edges is optional */
                if (edges != NULL) edges[k] = (edge){i, j};
                k++;
                /* actually a shifted linked list */
                if (!reachable(link, i, j)) swap(&link[i], &link[j]);
            }
        }
    }
    /*assert(k == nedges && "not enought edges CPLEX solution");*/
}

void assert_correctness(solution sol) {
    int visited = 0;
    int act = 0;
    do {
        visited++;
    } while ((act = sol->link[act]) != 0);

    assert(visited == sol->nedges);
}
