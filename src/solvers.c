#include "../include/solvers.h"

#include <assert.h>
#include <cplex.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

#include "../include/globals.h"
#include "../include/model_builder.h"
#include "../include/pqueue.h"
#include "../include/refinements.h"
#include "../include/tsp.h"
#include "../include/union_find.h"
#include "../include/utils.h"

void assert_correctness(solution sol);
solution TSPopt(instance inst, enum model_types model_type);
solution TSPminspantree(instance inst);
solution TSPgreedy(instance inst);
solution TSPgrasp(instance inst);
solution TSPextramileage(instance inst);
solution TSPvns(instance inst, int* succ);
solution TSPtabusearch(instance inst, int* succ);

solution solve(instance inst, enum model_types model_type) {
    /* initialize total wall-clock time of execution */
    struct timespec s, e;
    s.tv_sec = e.tv_sec = -1;
    stopwatch(&s, &e);

    solution sol;

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
            sol = TSPopt(inst, model_type);
            break;

        case MST:
            sol = TSPminspantree(inst);
            break;

        case GREEDY:
            sol = TSPgreedy(inst);
            break;

        case GRASP:
            sol = TSPgrasp(inst);
            break;

        case EXTRA_MILEAGE:
            sol = TSPextramileage(inst);
            break;

        case VNS:
            sol = TSPvns(inst, NULL);
            break;

        case TABU_SEACH:
            sol = TSPtabusearch(inst, NULL);
            break;

        case OPTIMAL_TOUR:
            sol = NULL;
            assert(model_type != OPTIMAL_TOUR &&
                   "tried to solve an optimal tour instance");
            break;
    }
    sol->solve_time = stopwatch(&s, &e);

    return sol;
}

solution TSPopt(instance inst, enum model_types model_type) {
    assert(inst != NULL);
    assert(inst->params != NULL && "no CPLEX params found");
    assert(inst->instance_type == TSP && "need TSP instance");

    /* create and populate solution */
    solution sol = create_solution(inst, model_type, inst->nnodes);
    sol->distance_time = compute_distmatrix(inst);

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
                           inst->params->timelimit * HF_INITIAL_PERC_TIME);
            /* alternative: work on branching node limit */
            /* CPXsetintparam(env, CPX_PARAM_NODELIM, 10); */

            /* recompute timelimit considering time spent on intial solution */
            inst->params->timelimit *= (1 - HF_INITIAL_PERC_TIME);

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
                           inst->params->timelimit * SF_INITIAL_PERC_TIME);
            /* alternative: work on branching node limit */
            CPXsetintparam(env, CPX_PARAM_NODELIM, 10);

            /* recompute timelimit considering time spent on intial solution */
            inst->params->timelimit *= (1 - SF_INITIAL_PERC_TIME);

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

    /* track deterministic time */
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

        case BENDERS: {
            struct timespec s, e;
            s.tv_sec = e.tv_sec = -1;
            stopwatch(&s, &e);
            perform_BENDERS(env, lp, inst, xstar, s, e);
            get_symmsol(xstar, sol);
            break;
        }

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
        case EXTRA_MILEAGE:
        case VNS:
        case TABU_SEACH:
            assert(0 == 1 && "tried to solve a metaheuristic");
            break;
    }

    free(xstar);

    /* retreive the time, in ticks and wall-clock */
    CPXgetdettime(env, &sol->end);

    /* retreive the min cost */
    CPXgetobjval(env, lp, &sol->zstar);
    /*sol->zstar = compute_zstar(inst, sol);*/

    /* add the solution to the pool associated with it's instance */
    add_solution(inst, sol);

    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);

    return sol;
}

solution TSPminspantree(instance inst) {
    assert(inst != NULL);
    assert(inst->instance_type == TSP && "need TSP instance");

    int nnodes = inst->nnodes;
    int nedges = nnodes * (nnodes - 1) / 2;

    solution sol = create_solution(inst, MST, nnodes);
    sol->distance_time = 0.0;

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

    /* refine it! */
    int* succ;
    if ((succ = edges_tosucc(sol->edges, nnodes)) == NULL) {
        print_error("solver didnt produced a tour");
    }
    sol->zstar += twoopt_refinement(inst, succ, nnodes);
    sol->zstar += threeopt_refinement(inst, succ, nnodes);
    /* store the succ as usual edges array */
    for (int j = 0; j < nnodes; j++) {
        sol->edges[j] = (edge){j, succ[j]};
    }
    free(succ);

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
        for (int j = 0; j < nnodes; j++) {
            if (visited[j]) continue;

            /* visit the node and store all the exiting edges */
            visited[j] = 1;
            for (int k = 0; k < uf->tns[j]->size; k++) {
                if (visited[uf->tns[j]->sons[k]]) continue;

                sol->edges[msti] = (edge){j, uf->tns[j]->sons[k]};

                msti++;
            }
        }
        free(visited);

        /* yeah i could reuse visitd and flip the bit, but this is executed just
         * one.. let me do this pls */
        int* edgecolors = (int*)calloc(sol->nedges, sizeof(int));
        for (int j = 1 + sol->nedges / 2; j < sol->nedges; j++) {
            edgecolors[j] = 1;
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

solution TSPgreedy(instance inst) {
    assert(inst != NULL);
    assert(inst->instance_type == TSP && "need TSP instance");

    int nnodes = inst->nnodes;

    /* track the best solution up to this point */
    solution sol = create_solution(inst, GREEDY, nnodes);
    sol->distance_time = 0.0;
    sol->zstar = DBL_MAX;

    /* succ used as visited: if -1 means not visited */
    int* succ = (int*)malloc(nnodes * sizeof(int));

    /* iterate over staring point (all possibilities) */
    for (int start = 0; start < 1; start++) {
        /* reset succ: -1 means not visited*/
        memset(succ, -1, nnodes * sizeof(int));

        if (VERBOSE) printf("[VERBOSE] greedy start %d\n", start + 1);

        /* add starting point to the actual tour */
        int act, next;
        act = start;
        double obj = 0.0;

        /* loop over nodes to visit: nodes -1 for start */
        int nunvisited = nnodes - 1;
        while (nunvisited--) {
            double weight = DBL_MAX;

            /* search for best new node */
            for (int i = 0; i < nnodes; i++) {
                if (succ[i] != -1) continue;
                if (i == act) continue;

                /* update best new edge */
                if (dist(act, i, inst) < weight) {
                    next = i;
                    weight = dist(act, i, inst);
                }
            }

            /* add the best new edge to the tour */
            succ[act] = next;
            obj += weight;

            if (EXTRA) printf("\tnext: %d, obj: %lf\n", next + 1, obj);

            act = next;
        }

        /* do not forget to close the loop! */
        succ[next] = start;
        obj += dist(next, start, inst);

        /* refine solution! (2opt only) */
        obj += twoopt_refinement(inst, succ, nnodes);

        if (VERBOSE) printf("\tfinish selecting, obj: %lf\n", obj);

        /* select best tour */
        if (obj < sol->zstar) {
            for (int i = 0; i < nnodes; i++) {
                sol->edges[i] = (edge){i, succ[i]};
            }
            sol->zstar = obj;
        }
    }
    /* refine it! (apply also 3opt) */
    sol->zstar += threeopt_refinement(inst, succ, nnodes);

    /* add the solution to the pool associated with it's instance */
    add_solution(inst, sol);

    free(succ);

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
    sol->distance_time = 0.0;
    sol->zstar = DBL_MAX;

    /* data structure for computation: visited flags, topk queue for performance
     * and tour container */
    topkqueue tk = topkqueue_create(GRASP_K);
    /* succ used as visited: if -1 means not visited */
    int* succ = (int*)malloc(nnodes * sizeof(int));

    /* iterate over staring point (randomly) until timelimit */
    double timeleft = inst->params->timelimit * (1.0 - GRASP_VNS_PERC_TIME);
    while (timeleft > 0) {
        /* initialize total wall-clock time for iteration */
        struct timespec s, e;
        s.tv_sec = e.tv_sec = -1;
        stopwatch_n(&s, &e);

        /* reset succ */
        memset(succ, -1, nnodes * sizeof(int));

        /* generate starting point to the actual tour */
        int start, act, next;
        start = act = rand_r(&seedp) % nnodes;
        double obj = 0.0;

        if (VERBOSE) printf("[VERBOSE] grasp start %d\n", start + 1);

        /* loop over nodes to visit: nodes -1 for start */
        int nunvisited = nnodes - 1;
        while (nunvisited--) {
            /* search for best new node */
            for (int i = 0; i < nnodes; i++) {
                if (succ[i] != -1) continue;
                if (i == act) continue; /* act is not visited yet */

                /* update top-k queue*/
                topkqueue_push(tk, dist(i, act, inst), i);
            }

            /* add the best new edge to the tour */
            next = topkqueue_randompick(tk);
            succ[act] = next;
            obj += dist(act, next, inst);

            if (EXTRA) printf("\tnext: %d, obj: %lf\n", next + 1, obj);

            act = next;
        }
        /* do not forget to close the loop! */
        succ[next] = start;
        obj += dist(next, start, inst);

        /* refine solution (2opt only)! */
        // obj += twoopt_refinement(inst, succ, nnodes);

        if (VERBOSE) printf(" -> obj: %lf\n", obj);

        /* select best tour */
        if (obj < sol->zstar) {
            for (int i = 0; i < nnodes; i++) {
                sol->edges[i] = (edge){i, succ[i]};
            }
            sol->zstar = obj;
        }

        /* update timelimit */
        timeleft -= stopwatch(&s, &e) / 1e3;
    }
    topkqueue_free(tk);

    /* refine it! */
    free(succ);
    if ((succ = edges_tosucc(sol->edges, nnodes)) == NULL) {
        print_error("solver didnt produced a tour");
    }
    free_solution(sol);
    /* TODO(lugot): AVOID so many frees */

    /* update timelimit and refine using vns */
    inst->params->timelimit *= GRASP_VNS_PERC_TIME;
    solution refined_sol = TSPvns(inst, succ);

    refined_sol->model_type = GRASP;
    free(succ);

    /* no need to add, added by TSPvns function */
    /* add_solution(inst, refined_sol); */

    return refined_sol;
}

solution TSPextramileage(instance inst) {
    assert(inst != NULL);
    assert(inst->instance_type == TSP && "need TSP instance");

    int nnodes = inst->nnodes;

    /* track the best solution up to this point */
    solution sol = create_solution(inst, EXTRA_MILEAGE, nnodes);
    sol->distance_time = 0.0;
    sol->zstar = 0.0;

    /* andew's monothone chain algorithm for convex */
    int k = 0;
    int* H = (int*)malloc(nnodes * 2 * sizeof(int));
    qsort(inst->nodes, nnodes, sizeof(struct node_t), nodelexcmp);

    for (int i = 0; i < nnodes; i++) {
        while (k >= 2) {
            node a, b, c;
            a = inst->nodes[H[k - 2]];
            b = inst->nodes[H[k - 1]];
            c = inst->nodes[i];
            if (ccw(a, b, c)) break;

            k--;
        }
        H[k++] = i;
    }
    for (int i = nnodes - 2, t = k + 1; i >= 0; i--) {
        while (k >= t) {
            node a, b, c;
            a = inst->nodes[H[k - 2]];
            b = inst->nodes[H[k - 1]];
            c = inst->nodes[i];
            if (ccw(a, b, c)) break;

            k--;
        }
        H[k++] = i;
    }
    int nhull = k; /* for EXTRA section */
    /* note: H closes the loop! */

    /* store the convex hull into the solution and track visited nodes */
    int* visited = (int*)calloc(nnodes, sizeof(int));
    for (int i = 0; i < k - 1; i++) {
        sol->edges[i] = (edge){H[i], H[i + 1]};
        sol->zstar += dist(H[i], H[i + 1], inst);

        visited[H[i]] = visited[H[i + 1]] = 1;
    }

    /* perform extra milage: select, among the remaining nodes, the one that
     * gives the lowest increment in the solution */
    k--; /* k conside the loop closed */
    int nunvisited = nnodes - k;
    while (nunvisited--) {
        double best_extra_milage = DBL_MAX;
        int next;  /* next unvisited node to pick */
        int edgei; /* edge index */

        /* iterate over unvisited nodes */
        for (int i = 0; i < nnodes; i++) {
            if (visited[i]) continue;

            /* iterate over saved edges */
            for (int j = 0; j < k; j++) {
                double extra_milage =
                    dist(i, sol->edges[j].i, inst) +
                    dist(i, sol->edges[j].j, inst) -
                    dist(sol->edges[j].i, sol->edges[j].j, inst);

                if (extra_milage < best_extra_milage) {
                    best_extra_milage = extra_milage;
                    next = i;
                    edgei = j;
                }
            }
        }

        /* add the new two edges to the saved ones:
         * actually perform a substitution and add the other one
         * also update the objective */
        int oldi = sol->edges[edgei].i;
        int oldj = sol->edges[edgei].j;
        sol->edges[edgei] = (edge){oldi, next};
        sol->edges[k] = (edge){oldj, next};
        k++;

        visited[next] = 1;

        sol->zstar += best_extra_milage;

        if (VERBOSE)
            printf("[VERBOSE] next: %d, break (%d, %d)\n", next + 1, oldi + 1,
                   oldj + 1);
    }

    /* refine it! */
    int* succ;
    if ((succ = edges_tosucc(sol->edges, nnodes)) == NULL) {
        print_error("solver didnt produced a tour");
    }
    sol->zstar += twoopt_refinement(inst, succ, nnodes);
    sol->zstar += threeopt_refinement(inst, succ, nnodes);
    /* store the succ as usual edges array */
    for (int i = 0; i < nnodes; i++) {
        sol->edges[i] = (edge){i, succ[i]};
    }
    free(succ);

    /* add the solution to the pool associated with it's instance */
    add_solution(inst, sol);

    if (EXTRA) {
        /* add the convex hull to the solution */
        sol->nedges += nhull - 1;
        sol->edges =
            (edge*)realloc(sol->edges, sol->nedges * sizeof(struct edge_t));

        int* edgecolors = (int*)calloc(sol->nedges, sizeof(int));
        int chi = sol->nedges - nhull + 1;
        for (int i = 0; i < nhull - 1; i++) {
            sol->edges[chi] = (edge){H[i], H[i + 1]};
            edgecolors[chi] = 1;

            chi++;
        }

        /* plot the instance with special code version 100 */
        plot_graphviz(sol, edgecolors, 100);
        free(edgecolors);

        /* return the correct solution */
        sol->nedges = sol->nedges - nhull + 1;
        sol->edges = realloc(sol->edges, sol->nedges * sizeof(struct edge_t));
    }

    free(visited);
    free(H);

    return sol;
}

solution TSPvns(instance inst, int* succ) {
    assert(inst != NULL);
    assert(inst->instance_type == TSP && "need TSP instance");

    int nnodes = inst->nnodes;
    srand(time(NULL));
    unsigned int seedp = time(NULL);

    /* track the best solution up to this point */
    solution sol = create_solution(inst, VNS, nnodes);
    sol->distance_time = 0.0;
    sol->zstar = DBL_MAX;

    /* if no successor array passed, create a random solution */
    int succ_tofree = 0;
    if (succ == NULL) {
        succ = randomtour(nnodes, seedp);
        succ_tofree = 1;
    }

    /* start the iteration! */
    int k = VNS_K_START;
    double timeleft = inst->params->timelimit;
    while (timeleft > 0 && k < VNS_K_MAX) {
        /* track the time for the current iteration */
        struct timespec s, e;
        s.tv_sec = e.tv_sec = -1;
        stopwatch(&s, &e);

        if (EXTRA) printf("[VERBOSE] kick size: %d\n", k);

        /* perturbe solution temp solution */
        kick(succ, nnodes, k);

        double obj = 0.0;
        for (int i = 0; i < nnodes; i++) obj += dist(i, succ[i], inst);

        if (EXTRA) printf("[VERBOSE] kicked objective: %lf\n", obj);

        /* find local optimum */
        printf("%lf\n", timeleft);
        obj += twoopt_refinement(inst, succ, nnodes);
        printf("%lf\n", timeleft);
        obj += threeopt_refinement(inst, succ, nnodes);
        printf("%lf\n", timeleft);

        if (EXTRA) printf("[VERBOSE] refined objective: %lf\n", obj);

        if (obj < sol->zstar) {
            /* store the succ as usual edges array */
            for (int i = 0; i < nnodes; i++) {
                sol->edges[i] = (edge){i, succ[i]};
            }

            if (EXTRA) {
                printf("\n[VERBOSE] Improved solution!\n");
                printf("[VEROBSE] -last opt: %lf, -new opt: %lf\n", sol->zstar,
                       obj);
                sleep(1);
            }
            /* update objective */
            sol->zstar = obj;

            /* restart k from initial value */
            k = VNS_K_MAX;

        } else {
            /* update k: increase kick strenght */
            k += VNS_K_STEP;
        }

        /* update timelimit */
        timeleft -= stopwatch(&s, &e) / 1000.0;
        printf("%lf\n", timeleft);
    }
    if (succ_tofree) free(succ);

    /* add the solution to the pool associated with it's instance */
    add_solution(inst, sol);

    return sol;
}

solution TSPtabusearch(instance inst, int* succ) {
    assert(inst != NULL);
    assert(inst->instance_type == TSP && "need TSP instance");

    int nnodes = inst->nnodes;
    srand(time(NULL));
    unsigned int seedp = time(NULL);

    /* track the best solution up to this point */
    solution sol = create_solution(inst, TABU_SEACH, nnodes);
    sol->distance_time = 0.0;
    sol->zstar = DBL_MAX;

    int* tabu_nodes = (int*)malloc(nnodes * sizeof(int));
    intset(tabu_nodes, -INF, nnodes);
    int tenure = TS_MAX_TENURE;

    /* if no successor array passed, create a random solution */
    int succ_tofree = 0;
    if (succ == NULL) {
        succ = randomtour(nnodes, seedp);
        succ_tofree = 1;
    }

    /* compute obj: it will be updated throught iterations */
    double obj = 0.0;
    for (int i = 0; i < nnodes; i++) obj += dist(i, succ[i], inst);

    /* start the iterations! */
    int k = 0; /* iteration counter */
    double timeleft = inst->params->timelimit;
    while (timeleft > 0) {
        /* track the time for the current iteration */
        struct timespec s, e;
        s.tv_sec = e.tv_sec = -1;
        stopwatch(&s, &e);

        int a, b;
        double delta =
            twoopt_tabu_pick(inst, succ, tabu_nodes, tenure, k, &a, &b);
        if (EXTRA) {
            printf("[EXTRA] iteration %d: delta %lf\n", k, delta);
        }

        /* actually perform the move, even if delta positive */
        twoopt_move(succ, nnodes, a, b);
        /* update the objective */
        obj += delta;

        if (delta > EPSILON) {
            /* climbing,  need to register nodes in tabulist */
            tabu_nodes[a] = tabu_nodes[b] = k;
        } else {
            /* downhill, update best solution if necessary */
            if (obj < sol->zstar) {
                /* store the succ as usual edges array */
                for (int i = 0; i < nnodes; i++) {
                    sol->edges[i] = (edge){i, succ[i]};
                }

                if (EXTRA) {
                    printf(
                        "[VERBOSE] Improved solution! last opt: %lf -> new "
                        "opt: %lf\n",
                        sol->zstar, obj);
                    // sleep(0.5);
                }
                /* update objective */
                sol->zstar = obj;
            }
        }

        /* update timelimit */
        timeleft -= stopwatch(&s, &e) / 1000.0;
        // printf("timeleft: %lf\n", timeleft);
        /* update iteration counter */
        k++;
    }
    if (succ_tofree) free(succ);

    /* add the solution to the pool associated with it's instance */
    add_solution(inst, sol);

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
            if (xstar[xxpos(i, j, sol->nedges)] > 0.5) {
                /* then store it in the solution*/
                sol->edges[k] = (edge){i, j};

                k++;
            }
        }
    }

    // TODO(any): check this
    /*assert(k == nedges && "not enought edges CPLEX solution");*/
}
