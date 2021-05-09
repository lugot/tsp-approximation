#include "../include/solvers.h"

#include <assert.h>
#include <cplex.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

#include "../include/adjlist.h"
#include "../include/globals.h"
#include "../include/model_builder.h"
#include "../include/pqueue.h"
#include "../include/tsp.h"
#include "../include/union_find.h"
#include "../include/utils.h"

void assert_correctness(solution sol);
solution TSPopt(instance inst, enum model_types model_type);
solution TSPminspantree(instance inst);
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

        case MST:
            return TSPminspantree(inst);
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

        case MST:
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
    sol->zstar = DBL_MAX;

    /* initialize total wall-clock time of execution and track deterministic
     * time */
    struct timespec s, e;
    s.tv_sec = e.tv_sec = -1;
    stopwatch(&s, &e);

    /* data structure for computation: visited flags, topk queue for performance
     * and tour container */
    int* visited = (int*)malloc(nnodes * sizeof(int));
    int* tour = (int*)malloc(nnodes * sizeof(int));

    /* iterate over staring point (all possibilities) */
    for (int start = 0; start < 1; start++) {
        /* reset visited */
        memset(visited, 0, nnodes * sizeof(int));

        if (VERBOSE) printf("[VERBOSE] greedy start %d\n", start);

        /* add starting point to the actual tour */
        int t = 0;
        tour[0] = start;
        visited[tour[0]] = 1;
        double obj = 0.0;

        t++;

        /* loop over nodes to visit */
        int next;
        int unvisited = nnodes - 1;
        while (unvisited > 0) {
            double weight = DBL_MAX;

            /* search for best new node */
            for (int i = 0; i < nnodes; i++) {
                if (visited[i]) continue;
                if (i == tour[t - 1]) continue;

                /* update best new edge */
                if (weight > dist(i, tour[t - 1], inst)) {
                    next = i;
                    weight = dist(i, tour[t - 1], inst);
                }
            }

            /* add the best new edge to the tour */
            tour[t] = next;
            obj += weight;
            visited[next] = 1;
            unvisited--;

            t++;

            if (EXTRA) printf("\tnext: %d, obj: %lf\n", next, obj);
        }

        /* do not forget to close the loop! */
        obj += dist(tour[t - 1], tour[1], inst);

        if (VERBOSE) printf("\tfinish selecting, obj: %lf\n", obj);

        /* select best tour */
        if (obj < sol->zstar) {
            for (int i = 0; i < nnodes; i++) {
                sol->edges[i] = (edge){tour[i], tour[(i + 1) % (nnodes)]};
            }
            sol->zstar = obj;
        }
    }
    /* track the solve time */
    sol->solve_time = stopwatch(&s, &e);

    /* add the solution to the pool associated with it's instance */
    add_solution(inst, sol);

    free(visited);
    free(tour);

    return sol;
}

solution TSPgrasp(instance inst) {
    assert(inst != NULL);
    assert(inst->instance_type == TSP && "need TSP instance");
    assert(inst->params != NULL);

    int nnodes = inst->nnodes;
    srand(time(NULL));
    unsigned int seedp = time(NULL);

    /* track the best solution up to this point */
    solution sol = create_solution(inst, GRASP, nnodes);
    sol->distance_time = compute_dist(inst);
    sol->zstar = DBL_MAX;

    /* data structure for computation: visited flags, topk queue for performance
     * and tour container */
    int* visited = (int*)malloc(nnodes * sizeof(int));
    topkqueue tk = topkqueue_create(GRASP_K);
    int* tour = (int*)malloc(nnodes * sizeof(int));

    /* iterate over staring point (randomly) until timelimit */
    while (inst->params->timelimit > 0) {
        /* initialize total wall-clock time for iteration */
        struct timespec s, e;
        s.tv_sec = e.tv_sec = -1;
        stopwatch_n(&s, &e);

        /* reset visited */
        memset(visited, 0, nnodes * sizeof(int));

        /* generate starting point to the actual tour */
        int t = 0;
        tour[0] = rand_r(&seedp) % nnodes;
        visited[tour[0]] = 1;
        double obj = 0.0;

        t++;

        if (VERBOSE) {
            printf("[VERBOSE] grasp start %d", tour[0]);
            fflush(0);
        }

        /* loop over nodes to visit */
        int unvisited = nnodes - 1;
        while (unvisited > 0) {
            /* search for best new node */
            for (int i = 0; i < nnodes; i++) {
                if (visited[i]) continue;
                if (i == tour[t - 1]) continue;

                /* update top-k queue*/
                topkqueue_push(tk, dist(i, tour[t - 1], inst), i);
            }

            /* add the best new edge to the tour */
            tour[t] = topkqueue_randompick(tk);
            obj += dist(tour[t - 1], tour[t], inst);
            visited[tour[t]] = 1;
            unvisited--;

            t++;

            if (EXTRA) printf("\tnext: %d, obj: %lf\n", tour[t - 1], obj);
        }

        /* do not forget to close the loop! */
        obj += dist(tour[t - 1], tour[0], inst);

        if (VERBOSE) printf(" -> obj: %lf\n", obj);

        /* select best tour */
        if (obj < sol->zstar) {
            for (int i = 0; i < nnodes; i++) {
                sol->edges[i] = (edge){tour[i], tour[(i + 1) % nnodes]};
            }
            sol->zstar = obj;
        }

        /* update timelimit */
        inst->params->timelimit -= stopwatch_n(&s, &e) / 1e9;
    }
    /* add the solution to the pool associated with it's instance */
    add_solution(inst, sol);

    topkqueue_free(tk);
    free(visited);
    free(tour);

    return sol;
}

solution TSPminspantree(instance inst) {
    assert(inst != NULL);
    assert(inst->instance_type == TSP && "need TSP instance");

    int nnodes = inst->nnodes;
    int nedges = nnodes * (nnodes - 1) / 2;

    solution sol = create_solution(inst, GREEDY, nnodes);
    sol->distance_time = compute_dist(inst);

    /* initialize total wall-clock time of execution time */
    struct timespec s, e;
    s.tv_sec = e.tv_sec = -1;
    stopwatch(&s, &e);

    /* build weighted edge list and sort it */
    wedge* wedges = (wedge*)malloc(nedges * sizeof(struct wedge_t));
    for (int i = 0; i < nnodes; i++) {
        for (int j = i + 1; j < nnodes; j++) {
            wedges[xpos(i, j, nnodes)] = (wedge){dist(i, j, inst), i, j};
        }
    }
    qsort(wedges, nedges, sizeof(struct wedge_t), wedgecmp);

    /* itearate over the edges and create the mst */
    union_find uf = uf_create(nnodes);
    for (int k = 0; k < nedges; k++) {
        int i, j;
        i = wedges[k].i;
        j = wedges[k].j;

        if (uf_same_set(uf, i, j)) continue;

        uf_union_set(uf, i, j);
    }

    /* save the solution: pre/post order visit the tree and save the nodes
     * encountered, starting from the uf root */
    sol->zstar = 0.0;
    int prev, next;
    prev = -1;
    int i = 0;
    while (uf_postorder(uf, &next)) {
        if (prev != -1) {
            sol->edges[i] = (edge){prev, next};
            sol->zstar += dist(prev, next, inst);

            i++;
        }

        prev = next;
    }
    /* do not forget to close the loop! */
    sol->edges[i] = (edge){prev, sol->edges[0].i};
    sol->zstar += dist(prev, sol->edges[0].i, inst);

    /* track the solve time */
    sol->solve_time = stopwatch(&s, &e);

    /* add the solution to the pool associated with it's instance */
    add_solution(inst, sol);

    if (EXTRA) {
        /* add the spanning tree to the solution */

        /* tree has nnodes -1 edges, realloc */
        sol->nedges += sol->nedges - 1;
        sol->edges =
            (edge*)realloc(sol->edges, sol->nedges * sizeof(struct edge_t));

        int msti = 1 + sol->nedges / 2; /* index of mst edge to add */
        int* visited = (int*)calloc(nnodes, sizeof(int));
        for (int i = 0; i < nnodes; i++) {
            if (visited[i]) continue;

            /* visit the node and store all the exiting edges */
            visited[i] = 1;
            for (int j = 0; j < uf->tns[i]->size; j++) {
                if (visited[uf->tns[i]->sons[j]]) continue;

                sol->edges[msti] = (edge){i, uf->tns[i]->sons[j]};

                msti++;
            }
        }
        free(visited);

        /* yeah i could reuse visitd and flip the bit, but this is executed just
         * one.. let me do this pls */
        int* edgecolors = (int*)calloc(sol->nedges, sizeof(int));
        for (int i = 1 + sol->nedges / 2; i < sol->nedges; i++) {
            edgecolors[i] = 1;
        }

        /* plot the instance with special code version 100 */
        plot_graphviz(sol, edgecolors, 100);
        free(edgecolors);

        /* return the correct solution */
        sol->nedges = 1 + sol->nedges / 2;
        sol->edges = realloc(sol->edges, sol->nedges * sizeof(struct edge_t));
    }

    free(wedges);
    uf_free(uf);

    return sol;
}

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
