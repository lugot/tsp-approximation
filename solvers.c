#define _GNU_SOURCE

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
solution TSPextramileage(instance inst);

double twoopt_refinement(instance inst, int* succ, int nnodes);
double twoopt_pick(instance inst, int* succ, int* a, int* b);
void twoopt_move(int* succ, int nnodes, int a, int b);
double threeopt_refinement(instance inst, int* succ, int nnodes);
double threeopt_pick(instance inst, int* succ, int* a, int* b, int* c);
void threeopt_move(int* succ, int nnodes, int a, int b, int c, instance inst);
solution VNS(instance inst, solution sol);

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

        case EXTRA_MILEAGE:
            return TSPextramileage(inst);
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
        case EXTRA_MILEAGE:
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

solution TSPminspantree(instance inst) {
    assert(inst != NULL);
    assert(inst->instance_type == TSP && "need TSP instance");

    int nnodes = inst->nnodes;
    int nedges = nnodes * (nnodes - 1) / 2;

    solution sol = create_solution(inst, MST, nnodes);
    sol->distance_time = 0.0;

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

    /* initialize total wall-clock time of execution and track deterministic
     * time */
    struct timespec s, e;
    s.tv_sec = e.tv_sec = -1;
    stopwatch(&s, &e);

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

        /* refine solution! */
        if (EXTRA) {
            /* save the solution before kick */
            for (int i = 0; i < nnodes; i++) {
                sol->edges[i] = (edge){i, succ[i]};
            }
            sol->inst = inst;
            /* plot the instance with special code version 100 */
            plot_graphviz(sol, NULL, 100);
        }
        kick(succ, nnodes, 5);

        if (VERBOSE) printf("\tfinish selecting, obj: %lf\n", obj);

        /* select best tour */
        if (obj < sol->zstar) {
            for (int i = 0; i < nnodes; i++) {
                sol->edges[i] = (edge){i, succ[i]};
            }
            sol->zstar = obj;
        }
    }
    /* track the solve time */
    sol->solve_time = stopwatch(&s, &e);

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
    while (inst->params->timelimit > 0) {
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

        /* refine solution! */
        /* twoopt_refinement(inst, succ, nnodes); */
        /* kick(succ, nnodes, 3); */

        if (VERBOSE) printf(" -> obj: %lf\n", obj);

        /* select best tour */
        if (obj < sol->zstar) {
            for (int i = 0; i < nnodes; i++) {
                sol->edges[i] = (edge){i, succ[i]};
            }
            sol->zstar = obj;
        }

        /* update timelimit */
        inst->params->timelimit -= stopwatch_n(&s, &e) / 1e9;
        break;
    }
    /* add the solution to the pool associated with it's instance */
    add_solution(inst, sol);

    topkqueue_free(tk);
    free(succ);

    return sol;
}

solution TSPextramileage(instance inst) {
    assert(inst != NULL);
    assert(inst->instance_type == TSP && "need TSP instance");

    int nnodes = inst->nnodes;

    /* track the best solution up to this point */
    solution sol = create_solution(inst, EXTRA_MILEAGE, nnodes);
    sol->distance_time = 0.0;
    sol->zstar = 0.0;

    /* initialize total wall-clock time of execution and track deterministic
     * time */
    struct timespec s, e;
    s.tv_sec = e.tv_sec = -1;
    stopwatch(&s, &e);

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
    while (nunvisited > 0) {
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
        nunvisited--;

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

    printf("Do VSN (zstar: %lf)\n",sol->zstar);
    sol = VNS(inst, sol);
    printf("Done (zstar: %lf)\n",sol->zstar);

    print_solution(sol, 1);

    /* track the solve time */
    sol->solve_time = stopwatch(&s, &e);

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

double twoopt_refinement(instance inst, int* succ, int nnodes) {
    double improvement = 0.0;

    int a, b;
    double delta = 1.0; /* using 0.0 as sentinel */
    a = b = 0;

    printf("succ: ");
    for (int i = 0; i < nnodes; i++) printf("%d ", succ[i]);
    printf("\n");

    /* iterate over 2opt moves until not improvable */
    while ((delta = twoopt_pick(inst, succ, &a, &b)) != 0.0) {
        if (EXTRA)
            printf("[VERBOSE] refinement on %d, %d delta %lf\n", a, b, delta);

        /* actually perform the move */
        twoopt_move(succ, nnodes, a, b);

        /* update the objective */
        improvement += delta; /* delta should be negative */
    }

    return improvement;
}

double twoopt_pick(instance inst, int* succ, int* a, int* b) {
    assert(inst != NULL);
    int nnodes, nedges;
    nnodes = nedges = inst->nnodes;

    double deltabest;
    deltabest = 0.0;
    *a = *b = 0;

    for (int i = 0; i < nnodes; i++) {
        for (int j = i + 1; j < nnodes; j++) {
            double delta = dist(i, j, inst) + dist(succ[i], succ[j], inst) -
                           (dist(i, succ[i], inst) + dist(j, succ[j], inst));

            if (delta < deltabest) {
                deltabest = delta;
                *a = i;
                *b = j;
            }
        }
    }

    return deltabest;
}

void twoopt_move(int* succ, int nnodes, int a, int b) {
    int aprime = succ[a], bprime = succ[b];

    /* b -> b' becomes b ~-> a' */
    reverse_path(succ, nnodes, aprime, b);
    /* a -> a' becomes a -> b */
    succ[a] = b;
    /* a' ~-> b becomes a' -> b' */
    succ[aprime] = bprime;
    /* b' ~-> a remains untouched */
    ;
}



double threeopt_refinement(instance inst, int* succ, int nnodes) {
    double improvement = 0.0;

    int a, b, c;
    double delta = 1.0; /* using 0.0 as sentinel */
    a = b = c = 0;

    /* iterate over 2opt moves until not improvable */
    while ((delta = threeopt_pick(inst, succ, &a, &b, &c)) != 0.0) {
        if (EXTRA) printf("[VERBOSE] refinement on %d, %d, %d delta %lf\n", a, b, c, delta);

        /* actually perform the move */
        threeopt_move(succ, nnodes, a, b, c, inst);

        /* update the objective */
        improvement += delta;  /* delta should be negative */
    }
    
    return improvement;
}

double threeopt_pick(instance inst, int* succ, int* a, int* b, int* c) {
    assert(inst != NULL);
    int nnodes, nedges;
    nnodes = nedges = inst->nnodes;

    double deltabest;
    deltabest = 0.0;
    *a = *b = *c = 0;

    for (int i = 0; i < nnodes; i++) {
        for (int j = i + 1; j < nnodes; j++) {
            for(int k = j + 1; k < nnodes; k++) {

                if(i == succ[j] || i == succ[k])
                    continue;
                if(j == succ[i] || j == succ[k])
                    continue;
                if(k == succ[i] || k == succ[j])
                    continue;
                
                double delta[4];
                for(int n=0; n<4; n++) delta[n] = - ( dist(i, succ[i], inst) + dist(j, succ[j], inst) + dist(k, succ[k], inst) );            

                int g = i;
                int swap = 0;
                int tj, tk;

                //makes k the third visited node
                while(succ[g]!=i){
                    if(swap == 0)
                    {
                        if(succ[g]==k) swap = 1;
                        if(succ[g]==j) swap = -1;
                    }
                    g = succ[g];
                }

                if(swap == 1)
                {
                    tk = j;
                    tj = k;
                }
                else
                {
                    tj = j;
                    tk = k;
                }
                
                //5 cases: preserving cases reducing to 2-opt and subtours
                delta[0] += dist(i, tj, inst) + dist(succ[i], tk, inst) +dist(succ[tj], succ[tk], inst); //i, j, succ[i], k, succ[j], succ[k]
                delta[1] += dist(i, succ[tj], inst) + dist(tk, succ[i], inst) + dist(tj, succ[tk], inst); //i, succ[j], k, succ[i], j, succ[k]
                delta[2] += dist(i, succ[tj], inst) + dist(tk, tj, inst) + dist(succ[i], succ[tk], inst); //i, succ[j], k, j, succ[i], succ[k]
                delta[3] += dist(i, tk, inst) + dist(succ[tj], succ[i], inst) + dist(tj, succ[tk], inst); //i, k, succ[j], succ[i], j, succ[k]


                for(int n=0; n<4; n++)
                    if (delta[n] < deltabest) {

                        deltabest = delta[n];
                        *a = i;
                        *b = tj;
                        *c = tk;
                    }
            }
        }
    }

    return deltabest;
}

void threeopt_move(int* succ, int nnodes, int a, int b, int c, instance inst) {
    int aprime = succ[a], bprime = succ[b], cprime = succ[c];

    double delta[4];

    int index = 0;

    delta[0] = dist(a, b, inst) + dist(succ[a], c, inst) +dist(succ[b], succ[c], inst); //i, j, succ[i], k, succ[j], succ[k]
    delta[1] = dist(a, succ[b], inst) + dist(c, succ[a], inst) + dist(b, succ[c], inst); //i, succ[j], k, succ[i], j, succ[k]
    delta[2] = dist(a, succ[b], inst) + dist(c, b, inst) + dist(succ[a], succ[c], inst); //i, succ[j], k, j, succ[i], succ[k]
    delta[3] = dist(a, c, inst) + dist(succ[b], succ[a], inst) + dist(b, succ[c], inst); //i, k, succ[j], succ[i], j, succ[k]

    double best = delta[0];

    for(int n=1; n<4; n++)
        if (delta[n] < best) {
            best = delta[n];
            index = n;
        }

    best = best - ( dist(a, succ[a], inst) + dist(b, succ[b], inst) + dist(c, succ[c], inst) );

    //if necessary swap tour from aprime to b
    if(index == 0 || index == 2){
        /* b -> b' becomes b ~-> a' */
        int* stack = (int*)malloc(nnodes * sizeof(int));
        int k = 0;

        stack[0] = aprime;
        k++;

        while (stack[k - 1] != b) {
            stack[k] = succ[stack[k - 1]];
            k++;
        }

        k--;
        while (k >= 1) {
            succ[stack[k]] = stack[k - 1];
            k--;
        }

        free(stack);
    }


    //if necessary swap tour from bprime to c
    if(index == 0 || index == 3){
        /* b -> b' becomes b ~-> a' */
        int* stack = (int*)malloc(nnodes * sizeof(int));
        int k = 0;

        stack[0] = bprime;
        k++;

        while (stack[k - 1] != c) {
            stack[k] = succ[stack[k - 1]];
            k++;
        }

        k--;
        while (k >= 1) {
            succ[stack[k]] = stack[k - 1];
            k--;
        }
        free(stack);
    }

    if(index == 0)
    {
        succ[a] = b;
        succ[aprime] = c;
        succ[bprime] = cprime;
    }
    if(index == 1)
    {
        succ[a] = bprime;
        succ[c] = aprime;
        succ[b] = cprime;
    }
    if(index == 2)
    {
        succ[a] = bprime;
        succ[c] = b;
        succ[aprime] = cprime;
    }
    if(index == 3)
    {
        succ[a] = c;
        succ[bprime] = aprime;
        succ[b] = cprime;
    }
}


void kick(int* succ, int nnodes, int strength) {
    assert(strength >= 2);
    assert(strength < nnodes / 2);

    /* tbh the worse line of code ive ever written */

    srand(time(NULL));

    /* phase 1: embedded successor array preparation */
    int embsize = 2 * strength; /* embedded size */

    int nselected = 0;
    int* selected = (int*)malloc(strength * sizeof(int));
    while (nselected < strength) {
        int candidate = rand() % nnodes;

        /* check if not previously selected or is at dist 1 from it */
        int valid = 1;
        for (int i = 0; i < nselected; i++) {
            if (candidate == selected[i] || candidate == succ[selected[i]] ||
                succ[candidate] == selected[i]) {
                valid = 0;
                break;
            }
        }

        if (valid) {
            selected[nselected] = candidate;
            nselected++;
        }
    }
    /* sort cause with reverses we're going to do a mess */
    int* pathlenghts = (int*)malloc(nnodes * sizeof(int));
    for (int i = 0; i < nselected; i++) {
        pathlenghts[selected[i]] = reachable(succ, 0, selected[i]);
    }
    qsort_r(selected, strength, sizeof(int), pathcmp, pathlenghts);

    /* with the nodes ordered, we can finally map nodes with successors to real
     * tour and delete the arcs from real tour */
    int* map = (int*)malloc(embsize * sizeof(int));
    for (int i = 0; i < nselected; i++) {
        /* even -> node, next odd -> its successor */
        map[2 * i] = selected[i];
        map[2 * i + 1] = succ[selected[i]];

        /* delete arcs in the real one */
        succ[selected[i]] = -1;
    }

    /* finally, create the embedding */
    int* embsucc = (int*)malloc(embsize * sizeof(int));
    for (int i = 0; i < embsize; i++) embsucc[i] = -1;

    /* to determine the target node while creating an arc, we have to be sure to
     * not create subtours (uf) and that the node is an actual target (relative
     * array) */
    union_find uf = uf_create(embsize);
    int* target = (int*)malloc(2 * strength * sizeof(int));
    for (int i = 0; i < embsize; i++) target[i] = 1; /* anyone is a target */

    /* set path sucecssor */
    for (int i = 0; i < 2 * strength; i++) {
        if (i % 2 != 1) continue; /* target odd position */

        embsucc[i] = (i + 1) % embsize;
        uf_union_set(uf, i, (i + 1) % embsize);
    }

    /* phase2: arc creation */
    int* buffer = (int*)malloc(embsize * sizeof(int));
    int tofix = strength;
    while (tofix--) { /* one arc creation for iteration */
        int a, b;
        int k;

        /* store in a buffer the legal starting nodes for a new arc */
        k = 0;
        for (int i = 0; i < embsize; i++) {
            if (embsucc[i] != -1) continue; /* legal == not a successor */

            buffer[k] = i;
            k++;
        }
        a = buffer[rand() % k];
        target[a] = 0; /* x cannot be targeted now */

        /* store in a buffer the legal ending position: we have to check that is
         * is a valid target and it's path is not the same as x */
        k = 0;
        for (int i = 0; i < embsize; i++) {
            /* in the last iteration there's clearly a single path */
            if (tofix != 0 && uf_same_set(uf, a, i)) continue;
            if (!target[i]) continue;

            buffer[k] = i;
            k++;
        }
        b = buffer[rand() % k];
        target[b] = 0;

        if (EXTRA) {
            printf("[VERBOSE] random selected %d(%d) -> %d(%d)\n", map[a] + 1,
                   a, map[a] + 1, b);
        }

        /* if b has not a successor, it means that we have to reverse the path
         * from a node u to b, but first we have to find this node u, the other
         * loose end of the path */

        if (embsucc[b] == -1) {
            int u;

            /* pick the fartherst node from b */
            int maxpathlen = 0;
            for (int i = 0; i < embsize; i++) {
                if (!uf_same_set(uf, b, i)) continue;

                int pathlen =
                    max(reachable(embsucc, b, i), reachable(embsucc, i, b));

                if (pathlen > maxpathlen) {
                    maxpathlen = pathlen;
                    u = i;
                }
            }

            if (VERBOSE)
                printf("[VERBOSE] perform reverse %d(%d) ~-> %d(%d)\n",
                       map[u] + 1, u, map[b] + 1, b);

            /* u becomes a loose end (in opposite direction) */
            reverse_path(embsucc, embsize, u, b);
            embsucc[u] = -1;
            target[u] = 1;

            reverse_path(succ, nnodes, map[u], map[b]);
            succ[map[u]] = -1;
        }

        /* create the new arc: set both embedded and real successor array and
         * update the union find data structure */
        embsucc[a] = b;
        uf_union_set(uf, a, b);
        succ[map[a]] = map[b];
    }

    /* *looking at the heap an being vito corleone* look how they massacred my
     * boy*/
    uf_free(uf);
    free(pathlenghts);
    free(selected);
    free(buffer);
    free(embsucc);
    free(map);
}

solution VNS(instance inst, solution sol)
{
	int nnodes = inst->nnodes;
	solution tmp_sol = create_solution(inst, EXTRA_MILEAGE, nnodes);

	for (int i = 0; i < nnodes; i++) {
    	tmp_sol->edges[i] = sol->edges[i];
    }
    	tmp_sol->zstar = sol->zstar ;

    /* refine it! */
    int* succ;
   	if ((succ = edges_tosucc(tmp_sol->edges, nnodes)) == NULL) {
   	    print_error("solver didnt produced a tour");
    }

    int starting_value = 5;
    int kmax = 20;

    int k = starting_value;

    while (inst->params->timelimit > 0 && k<kmax) {
        /* track the time for the current iteration */
        struct timespec s, e;
        s.tv_sec = e.tv_sec = -1;
        stopwatch(&s, &e);

        if(EXTRA) 
        	printf("\n[VERBOSE] Kick size: %d\n", k);
        /* perturbe solution temp solution */
        kick(succ, nnodes, k);

        double obj = 0;

        for(int i=0; i<nnodes; i++) obj+=dist(i, succ[i], inst);

        if(EXTRA) 
        	printf("\n[VERBOSE] Kicked Objective: %lf\n",obj);

        /* find local optimum */
    	obj += twoopt_refinement(inst, succ, nnodes);
    	obj += threeopt_refinement(inst, succ, nnodes);

    	if(EXTRA) 
    		printf("\n[VERBOSE] Refined Objective: %lf\n",obj);

    	tmp_sol->zstar = obj;
    	/* store the succ as usual edges array */
    	for (int i = 0; i < nnodes; i++) {
    	    tmp_sol->edges[i] = (edge){i, succ[i]};

    	    /* if new global optimum save solution */
    	    if(obj<sol->zstar)
    	    	sol->edges[i] = (edge){i, succ[i]};    	    	
    	}

    	if(obj<sol->zstar){

    	    if(EXTRA) {
    	    	printf("\n[VERBOSE] Improved solution!\n");
    	    	printf("[VEROBSE] -last opt: %lf, -new opt: %lf\n", sol->zstar, obj);
    	    	sleep(5);
    	    }

    		sol->zstar = obj;
    	   	k = starting_value;
    	   	
    	}
    	else k++;

        /* update timelimit */
        inst->params->timelimit -= stopwatch(&s, &e) / 1000.0;
    }

    free(succ);
    return sol;
}