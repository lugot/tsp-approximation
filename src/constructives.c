#include "../include/constructives.h"

#include <assert.h>
#include <float.h>
#include <string.h>
#include <time.h>

#include "../include/globals.h"
#include "../include/pqueue.h"
#include "../include/utils.h"

solution TSPgreedy(instance inst) {
    assert(inst != NULL);

    int nnodes = inst->nnodes;

    /* track the best solution up to this point */
    solution sol = create_solution(inst, GREEDY, nnodes);
    sol->distance_time = 0.0;
    sol->zstar = DBL_MAX;

    /* succ used as visited: if -1 means not visited */
    int* succ = (int*)malloc(nnodes * sizeof(int));

    /* initialize total wall-clock time */
    struct timespec s, e;
    s.tv_sec = e.tv_sec = -1;
    stopwatch(&s, &e);

    /* iterate over staring point (all possibilities) */
    for (int start = 0;
         start < nnodes && stopwatch(&s, &e) / 1000.0 < inst->params->timelimit;
         start++) {
        /* reset succ: -1 means not visited*/
        memset(succ, -1, nnodes * sizeof(int));

        if (EXTRA_VERBOSE) printf("[VERBOSE] greedy start %d\n", start + 1);

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

            if (EXTRA_VERBOSE) printf("\tnext: %d, obj: %lf\n", next + 1, obj);

            act = next;
        }

        /* do not forget to close the loop! */
        succ[next] = start;
        obj += dist(next, start, inst);

        if (EXTRA_VERBOSE) printf("\tfinish selecting, obj: %lf\n", obj);

        /* select best tour */
        if (obj < sol->zstar) {
            for (int i = 0; i < nnodes; i++) {
                sol->edges[i] = (edge){i, succ[i]};
            }
            sol->zstar = obj;

            tracker_add(sol->t, stopwatch(&s, &e), obj);
        }
    }

    free(succ);

    return sol;
}

solution TSPgrasp(instance inst, int onesolution) {
    assert(inst != NULL);
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

    /* initialize total wall-clock time */
    struct timespec s, e;
    s.tv_sec = e.tv_sec = -1;
    stopwatch(&s, &e);

    /* iterate over staring point (randomly) until timelimit */
    while (stopwatch(&s, &e) / 1000.0 < inst->params->timelimit) {
        /* reset succ */
        memset(succ, -1, nnodes * sizeof(int));

        /* generate starting point to the actual tour */
        int start, act, next;
        start = act = rand_r(&seedp) % nnodes;
        double obj = 0.0;

        if (EXTRA_VERBOSE) printf("[VERBOSE] grasp start %d\n", start + 1);

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

            if (EXTRA_VERBOSE) printf("\tnext: %d, obj: %lf\n", next + 1, obj);

            act = next;
        }
        /* do not forget to close the loop! */
        succ[next] = start;
        obj += dist(next, start, inst);

        if (EXTRA_VERBOSE) printf(" -> obj: %lf\n", obj);

        /* select best tour */
        if (obj < sol->zstar) {
            for (int i = 0; i < nnodes; i++) {
                sol->edges[i] = (edge){i, succ[i]};
            }
            sol->zstar = obj;

            tracker_add(sol->t, stopwatch(&s, &e), obj);
        }

        if (onesolution) break;
    }

    topkqueue_free(tk);
    free(succ);

    return sol;
}

solution TSPextramileage(instance inst) {
    assert(inst != NULL);

    int nnodes = inst->nnodes;

    /* track the best solution up to this point */
    solution sol = create_solution(inst, EXTRA_MILEAGE, nnodes);
    sol->distance_time = 0.0;
    sol->zstar = 0.0;

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
    while (nunvisited--) {
        double best_extra_milage = DBL_MAX;
        int next;  /* next unvisited node to pick */
        int edgei; /* edge index */

        int kk = 0;

        /* iterate over unvisited nodes */
        for (int i = 0; i < nnodes; i++) {
            if (visited[i]) continue;

            /* iterate over saved edges */
            for (int j = 0; j < k; j++) {
                kk++;
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

    if (EXTRA_VERBOSE) {
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
