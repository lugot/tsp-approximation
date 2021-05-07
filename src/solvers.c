#include "../include/solvers.h"

#include <assert.h>
#include <cplex.h>
#include <float.h>
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
solution TSPopt(instance inst, enum model_types model_type);
solution TSPgreedy(instance inst);
solution TSPgrasp(instance inst);

solution solve(instance inst, enum model_types model_type) {
    switch (model_type) {
        case NOSEC:
        case MTZ_STATIC:
        case MTZ_LAZY:
        case GGLIT_STATIC:
        case GGLECT_STATIC:
        case GGLIT_LAZY:
        case BENDERS:
        case BENDERS_CALLBACK:
        case HARD_FIXING:
        case SOFT_FIXING:
            return TSPopt(inst, model_type);
            break;

        case GREEDY:
            return TSPgreedy(inst);
            break;

        case GRASP:
            return TSPgrasp(inst);
            break;

        case OPTIMAL_TOUR:
            assert(model_type != OPTIMAL_TOUR &&
                   "tried to solve an optimal tour instance");
            break;
    }
    return NULL;
}

solution TSPopt(instance inst, enum model_types model_type) {
    assert(inst != NULL);
    assert(inst->params != NULL && "no CPLEX params found");
    assert(inst->instance_type == TSP && "need TSP instance");

    /* create and populate solution */
    solution sol = create_solution(inst, model_type, inst->nnodes);
    sol->distance_time = compute_dist(inst);

    /* open CPLEX model */
    int error;
    CPXENVptr env = CPXopenCPLEX(&error);
    CPXLPptr lp = CPXcreateprob(env, &error, "TSP");

    /* set params:
     * - timelimit (timetype = 0 -> seconds, timetype = 1 -> ticks)
     * - EpInteger: CPLEX tollerance to declare a variable integer, important in
     * bigM
     * - EpRightHandSide: the less or equal satisfied up to this tollerance,
     * usually 1e-5, with bigM 1e-9 */
    CPXsetdblparam(env, CPX_PARAM_TILIM, inst->params->timelimit);
    CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0);
    CPXsetdblparam(env, CPX_PARAM_EPRHS, 1e-9);
    // CPXsetintparam(env, CPX_PARAM_RANDOMSEED, seed); todo: add a seed to
    // Cplexrun (Passed by TSP opt?) TODO

    /* save the log of the execution */
    char logfile[] = "execution.log";
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
            /* set iniitial fraction of timelimit for first execution */
            CPXsetdblparam(env, CPX_PARAM_TILIM,
                           inst->params->timelimit / HF_INITIAL_FRACTION_TIME);
            /* alternative: work on branching node limit */
            /* CPXsetintparam(env, CPX_PARAM_NODELIM, 10); */

            /* recompute timelimit considering time spent on intial solution */
            inst->params->timelimit *= (1 - 1.0 / HF_INITIAL_FRACTION_TIME);

            /* set benders callbacks as blackbox for matheuristic */
            CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE;
            if (CPXcallbacksetfunc(env, lp, contextid,
                                   add_BENDERS_sec_callback_driver, sol)) {
                print_error("CPXcallbacksetfunc() error");
            }
            break;
        }

        case SOFT_FIXING: {
            /* set iniitial fraction of timelimit for first execution */
            CPXsetdblparam(env, CPX_PARAM_TILIM,
                           inst->params->timelimit / SF_INITIAL_FRACTION_TIME);
            /* alternative: work on branching node limit */
            CPXsetintparam(env, CPX_PARAM_NODELIM, 10);

            /* recompute timelimit considering time spent on intial solution */
            inst->params->timelimit *= (1 - 1.0 / SF_INITIAL_FRACTION_TIME);

            /* set benders callbacks as blackbox for matheuristic */
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

    /* initialize total wall-clock time of execution and track deterministic
     * time */
    struct timespec s, e;
    s.tv_sec = e.tv_sec = -1;
    stopwatch(&s, &e);
    CPXgetdettime(env, &sol->start);

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
            perform_BENDERS(env, lp, inst, xstar, s, e);
            get_symmsol(xstar, sol);
            break;

        case HARD_FIXING:
            perform_HARD_FIXING(env, lp, inst, xstar, HF_PERCENTAGE);
            get_symmsol(xstar, sol);
            break;

        case SOFT_FIXING:
            perform_SOFT_FIXING(env, lp, inst, xstar);
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

        case GREEDY:
        case GRASP:
            assert(model_type != GREEDY && model_type != GRASP &&
                   "tried to solve a metaheuristic");
            break;
    }

    free(xstar);

    /* retreive the time, in ticks and wall-clock */
    CPXgetdettime(env, &sol->end);
    sol->solve_time = stopwatch(&s, &e);

    /* retreive the min cost */
    CPXgetobjval(env, lp, &sol->zstar);
    /*sol->zstar = compute_zstar(inst, sol);*/

    /* add the solution to the pool associated with it's instance */
    add_solution(inst, sol);

    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);

    return sol;
}

solution TSPgreedy(instance inst) {
    assert(inst != NULL);
    assert(inst->instance_type == TSP && "need TSP instance");

    int nnodes = inst->nnodes;

    /* track the best solution up to this point */
    solution sol = create_solution(inst, GREEDY, nnodes);
    sol->distance_time = compute_dist(inst);
    print_instance(inst, 1);

    /* initialize total wall-clock time of execution and track deterministic
     * time */
    struct timespec s, e;
    s.tv_sec = e.tv_sec = -1;
    stopwatch(&s, &e);

    /* track visited nodes for each starting point (will be memsetted) */
    int* visited = (int*)malloc(nnodes * sizeof(int));

    /* track best solution so far */
    int* best_tour = (int*)malloc((nnodes - 1) * sizeof(int));
    double best_obj = DBL_MAX;

    /* iterate over staring point (all possibilities) */
    int* tour = (int*)malloc((nnodes - 1) * sizeof(int));
    double obj;
    for (int start = 0; start < nnodes; start++) {
        /* reset visited */
        memset(visited, 0, nnodes * sizeof(int));

        if (VERBOSE) printf("[VERBOSE] greedy start %d\n", start);

        /* add starting point to the actual tour */
        int t = 0;
        tour[t] = start;
        t++;

        visited[start] = 1;
        obj = 0.0;

        /* loop over nodes to visit */
        int unvisited = nnodes - 1;
        int next;
        while (unvisited > 0) {
            double weight = DBL_MAX; /* hopefully cached */

            /* search for best new node */
            for (int i = 0; i < nnodes; i++) {
                if (visited[i]) continue;
                if (i == tour[t - 1]) continue;

                /* update best new edge */
                if (weight > dist(i, tour[t-1], inst)) {
                    next = i;
                    weight = dist(i, tour[t-1], inst);
                }
            }

            /* add the best new edge to the tour */
            obj += weight;
            tour[t] = next;
            t++;

            visited[next] = 1;
            unvisited--;

            if (VERBOSE) printf("\tnext: %d, obj: %lf\n", next, obj);
        }

        /* do not forget to close the loop! */
        obj += inst->adjmatrix[next][start];

        if (VERBOSE) printf("\tfinish selecting, obj: %lf\n", obj);

        /* select best tour */
        if (obj < best_obj) {
            memcpy(best_tour, tour, (nnodes - 1) * sizeof(int));
            best_obj = obj;
        }
    }

    /* store the best tour in the solution */
    for (int i = 0; i < nnodes; i++) {
        sol->edges[i] = (edge){tour[i], tour[(i+1) % nnodes]};
    }
    sol->zstar = best_obj;

    /* track the solve time */
    sol->solve_time = stopwatch(&s, &e);

    /* add the solution to the pool associated with it's instance */
    add_solution(inst, sol);

    free(visited);
    free(tour);
    free(best_tour);

    return sol;
}

solution TSPgrasp(instance inst) { return NULL; }

void get_symmsol(double* xstar, solution sol) {
    /* index over selected edges */
    int k = 0;

    /* check for edges j>i if is selected */
    for (int i = 0; i < sol->nedges; i++) {
        for (int j = i + 1; j < sol->nedges; j++) {
            if (xstar[xpos(i, j, sol->nedges)] > 0.5) {
                /* then store it in the solution*/
                sol->edges[k] = (edge){i, j};

                k++;
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
                sol->edges[k] = (edge){i, j};

                k++;
            }
        }
    }

    // TODO(any): check this
    /*assert(k == nedges && "not enought edges CPLEX solution");*/
}
