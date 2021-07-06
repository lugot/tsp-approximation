#include "../include/metaheuristics.h"

#include <assert.h>
#include <float.h>
#include <time.h>
#include <unistd.h>

#include "../include/globals.h"
#include "../include/refinements.h"
#include "../include/utils.h"

solution TSPvns(instance inst, int* succ) {
    assert(inst != NULL);

    int nnodes = inst->nnodes;
    srand(time(NULL));
    unsigned int seedp = time(NULL);

    /* track the best solution up to this point */
    solution sol = create_solution(inst, VNS_RANDOM, nnodes);
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

        if (EXTRA_VERBOSE) printf("[VERBOSE] kick size: %d\n", k);

        /* perturbe solution temp solution */
        kick(succ, nnodes, k);

        double obj = 0.0;
        for (int i = 0; i < nnodes; i++) obj += dist(i, succ[i], inst);

        if (EXTRA_VERBOSE) printf("[VERBOSE] kicked objective: %lf\n", obj);

        /* find local optimum */
        printf("%lf\n", timeleft);
        obj += twoopt_refinement(inst, succ, nnodes);
        printf("%lf\n", timeleft);
        obj += threeopt_refinement(inst, succ, nnodes);
        printf("%lf\n", timeleft);

        if (EXTRA_VERBOSE) printf("[VERBOSE] refined objective: %lf\n", obj);

        if (obj < sol->zstar) {
            /* store the succ as usual edges array */
            for (int i = 0; i < nnodes; i++) {
                sol->edges[i] = (edge){i, succ[i]};
            }

            if (EXTRA_VERBOSE) {
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

    return sol;
}

solution TSPtabusearch(instance inst, int* succ) {
    assert(inst != NULL);

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
        if (EXTRA_VERBOSE) {
            printf("[VERBOSE] iteration %d: delta %lf\n", k, delta);
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

                if (EXTRA_VERBOSE) {
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

    return sol;
}
