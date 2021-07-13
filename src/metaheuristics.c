#include "../include/metaheuristics.h"

#include <assert.h>
#include <float.h>
#include <time.h>
#include <unistd.h>

#include "../include/globals.h"
#include "../include/refinements.h"
#include "../include/utils.h"

solution child(instance inst, solution parent1, solution parent2);
solution TSPrandom(instance inst);

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

    /* initialize total wall-clock time */
    struct timespec s, e;
    s.tv_sec = e.tv_sec = -1;
    stopwatch(&s, &e);

    /* start the iteration! */
    int k = VNS_K_START;
    int first_iter = 1;
    while (stopwatch(&s, &e) / 1000.0 < inst->params->timelimit &&
           k < VNS_K_MAX) {
        if (EXTRA_VERBOSE) printf("[VERBOSE] kick size: %d\n", k);

        /* perturbe solution temp solution */
        if (!first_iter)
            kick(succ, nnodes, k);
        else
            first_iter = 0;

        double obj = 0.0;
        for (int i = 0; i < nnodes; i++) obj += dist(i, succ[i], inst);

        if (EXTRA_VERBOSE) printf("[VERBOSE] kicked objective: %lf\n", obj);

        /* find local optimum */
        obj += twoopt_refinement(inst, succ, nnodes, &s, &e);
        /* obj += threeopt_refinement(inst, succ, nnodes); */

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
            }
            /* update objective */
            sol->zstar = obj;

            /* restart k from initial value */
            k = VNS_K_MAX;

        } else {
            /* update k: increase kick strenght */
            k += VNS_K_STEP;
        }
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
    solution sol = create_solution(inst, TABU_SEACH_RANDOMSTART, nnodes);
    sol->distance_time = 0.0;
    sol->zstar = DBL_MAX;

    int* tabu_nodes = (int*)malloc(nnodes * sizeof(int));
    intset(tabu_nodes, -INF, nnodes);
    int tenure = TS_MAX_TENURE;
    int diversification = 1;
    int downhill = 1;

    /* if no successor array passed, create a random solution */
    int succ_tofree = 0;
    if (succ == NULL) {
        succ = randomtour(nnodes, seedp);
        succ_tofree = 1;
    }

    /* initialize total wall-clock time */
    struct timespec s, e;
    s.tv_sec = e.tv_sec = -1;
    stopwatch(&s, &e);

    /* compute obj: it will be updated throught iterations */
    double obj = 0.0;
    for (int i = 0; i < nnodes; i++) obj += dist(i, succ[i], inst);

    /* start the iterations! */
    int k = 0; /* iteration counter */
    while (stopwatch(&s, &e) / 1000.0 < inst->params->timelimit) {
        int a, b;
        double delta =
            twoopt_tabu_pick(inst, succ, tabu_nodes, tenure, k, &a, &b);
        if (EXTRA_VERBOSE) {
            printf("[VERBOSE] iteration %d: delta %lf\n", k, delta);
        }

        if (delta > EPSILON && downhill) {
            /* local optimum, save */
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
                }
                /* update objective */
                sol->zstar = obj;

                tracker_add(sol->t, stopwatch(&s, &e), sol->zstar);
            }

            downhill = 0;
        }
        if (delta < -EPSILON) downhill = 1;

        /* actually perform the move, even if delta positive */
        twoopt_move(succ, nnodes, a, b);
        /* update the objective */
        obj += delta;

        if (delta > EPSILON) {
            /* climbing,  need to register nodes in tabulist */
            tabu_nodes[a] = tabu_nodes[b] = k;
        }
        /* else downhill, but do nothing, save sol only in local optimum */

        /* update iteration counter */
        k++;

        if ((k + 1) % TS_PHASEDURATION == 0) {
            diversification = (diversification + 1) % 2;

            if (diversification)
                tenure = TS_MIN_TENURE;
            else
                tenure = TS_MAX_TENURE;

            if (VERBOSE) {
                printf(
                    "[VERBOSE] iteration %d, diversification = %d, tenure %d\n",
                    k, diversification, tenure);
            }
        }
    }

    /* maybe last donwhill was the best one, save here! */
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
        }
        /* update objective */
        sol->zstar = obj;

        tracker_add(sol->t, stopwatch(&s, &e), sol->zstar);
    }

    if (succ_tofree) free(succ);
    free(tabu_nodes);

    return sol;
}

solution TSPgenetic(instance inst) {
    assert(inst != NULL);
    if (EXTRA_VERBOSE) {
        printf("[VERBOSE] Solving using genetic algorithm!\n");
    }

    /* initialize total wall-clock time of execution time */
    struct timespec s, e;
    s.tv_sec = e.tv_sec = -1;
    stopwatch(&s, &e);

    int nnodes = inst->nnodes;
    solution best = create_solution(inst, GENETIC, nnodes);

    double timeleft = inst->params->timelimit;

    int n, k;
    if (nnodes < 100) {
        n = 1000;
        k = 100;
    } else if (nnodes < 400) {
        n = 600;
        k = 60;
    } else if (nnodes < 750) {
        n = 120;
        k = 40;
    } else {
        /* no refinement for larger instances */
        n = 1000;
        k = 100;
        timeleft *= 0.9;
        /* time for final refinement on very big instance */
    }

    solution* population = (solution*)calloc(n + k, sizeof(solution));
    int* alive = (int*)calloc(n + k, sizeof(int));

    /* track the time for the current iteration */
    struct timespec s_left, e_left;
    s_left.tv_sec = e_left.tv_sec = -1;
    stopwatch(&s_left, &e_left);

    for (int i = 0; i < n; i++) {
        population[i] = TSPrandom(inst);
        alive[i] = 1;
    }

    /* update timelimit */
    timeleft -= stopwatch(&s_left, &e_left) / 1000.0;

    if (EXTRA_VERBOSE) printf("[VERBOSE] population created!\n");

    while (timeleft > 0) {
        s_left.tv_sec = e_left.tv_sec = -1;
        stopwatch(&s_left, &e_left);

        // not needed: All become 1 at the
        // end of the children creation
        int* indexes = (int*)calloc(k, sizeof(int));

        int index = 0;
        int current = 0;

        // create k children
        while (current < n + k) {
            if (alive[current] == 0) {
                // select randomly 2 parents
                solution parent1, parent2;
                int done = 0;

                int parent1_idx;
                while (!done) {
                    int pick = rand() % (n + k);
                    if (pick == current) continue;
                    if (alive[pick] == 1) {
                        parent1 = population[pick];
                        parent1_idx = pick;
                        done = 1;
                    }
                }

                done = 0;
                while (!done) {
                    int pick = rand() % (n + k);
                    if (pick == current) continue;
                    if (pick == parent1_idx) continue;
                    /* non oso immaginare i casini in caso di collisione */

                    if (alive[pick] == 1) {
                        parent2 = population[pick];
                        done = 1;
                    }
                }

                indexes[index] = current;
                /* free_solution(population[current]); */
                population[current] = child(inst, parent1, parent2);

                index++;
            }

            current++;
        }

        for (int j = 0; j < k; j++) {
            alive[indexes[j]] = 1;
        }

        double min = DBL_MAX;
        int min_idx;
        for (int x = 0; x < n + k; x++) {
            if (population[x]->zstar < min) {
                if (EXTRA_VERBOSE) printf("Best solution: %lf\n", min);
                min = population[x]->zstar;
                min_idx = x;
            }
        }

        tracker_add(best->t, stopwatch(&s, &e), population[min_idx]->zstar);

        // keep best n solutions
        for (int j = 0; j < k; j++) {
            double max = -1;
            int pick = -1;

            for (int x = 0; x < n + k; x++) {
                if (population[x]->zstar > max) {
                    max = population[x]->zstar;
                    pick = x;
                }
            }

            alive[pick] = 0;
        }

        /* update timelimit */
        timeleft -= stopwatch(&s_left, &e_left) / 1000.0;
        free(indexes);
    }

    // return best solution
    double min = DBL_MAX;
    int pick = -1;

    for (int i = 0; i < n + k; i++) {
        if (population[i]->zstar < min) {
            min = population[i]->zstar;
            pick = i;
        }
    }

    best = population[pick];

    /*if never refined, refine it! */
    if (nnodes >= 750) {
        int* succ = (int*)calloc(nnodes, sizeof(int));
        if ((succ = edges_tosucc(best->edges, nnodes)) == NULL) {
            print_error("solver didnt produced a tour");
        }

        s.tv_sec = e.tv_sec = -1;
        stopwatch(&s, &e);
        inst->params->timelimit *= 0.1;

        best->zstar += twoopt_refinement(inst, succ, nnodes, &s, &e);
        inst->params->timelimit /= 0.1;

        /* store the succ as usual edges array */
        for (int i = 0; i < nnodes; i++) {
            best->edges[i] = (edge){i, succ[i]};
        }
        free(succ);
    }

    /* track the solve time */
    best->solve_time = stopwatch(&s, &e);

    /* add the solution to the pool associated with it's instance */
    add_solution(inst, best);

    free(alive);
    for (int i = 0; i < n + k; i++) {
        if (i == pick) continue;
        /* free_solution(population[i]); */
    }
    free(population);

    return best;
}

solution child(instance inst, solution parent1, solution parent2) {
    if (EXTRA_VERBOSE) {
        printf("[VERBOSE] creating children!\n");
    }

    int nnodes = inst->nnodes;
    int split = nnodes / 2;

    solution child = create_solution(inst, GENETIC, inst->nnodes);

    int* visited = (int*)calloc(nnodes, sizeof(int));

    int* chromosome = (int*)calloc(nnodes, sizeof(int));
    int* succ = (int*)calloc(nnodes, sizeof(int));
    succ = edges_tosucc(parent1->edges, nnodes);

    /* nodes included in the solution */
    int visnodes = 0;

    chromosome[0] = 0;
    visited[0] = 1;
    visnodes++;

    for (int i = 1; i < split; i++) {
        chromosome[i] = succ[chromosome[i - 1]];
        visited[chromosome[i]] = 1;
        visnodes++;
    }

    succ = edges_tosucc(parent2->edges, nnodes);

    int done = 0;
    // aux skip visited nodes
    int aux = 0;

    while (!done) {
        if (visited[split + aux] == 1) {
            aux++;
        } else {
            chromosome[split] = split + aux;
            visited[chromosome[split]] = 1;
            aux = 0;
            visnodes++;
            done = 1;
        }
    }

    aux = 0;

    for (int i = split + 1; i < nnodes && aux < nnodes - split; i++) {
        int next = succ[chromosome[i - 1]];

        while (visited[next] == 1) {
            next = succ[next];
            aux++;
        }

        chromosome[i] = next;
        visited[next] = 1;
        visnodes++;
    }

    // chromosome contains a <nnodes cycle to be completed with extr mileage
    for (int i = 0; i < visnodes; i++) {
        child->edges[i] = (edge){chromosome[i], chromosome[(i + 1) % visnodes]};
        child->zstar +=
            dist(chromosome[i], chromosome[(i + 1) % visnodes], inst);
    }

    /* extramileage */
    int nunvisited = nnodes - visnodes;
    while (nunvisited--) {
        double best_extra_milage = DBL_MAX;
        int next;  /* next unvisited node to pick */
        int edgei; /* edge index */

        /* iterate over unvisited nodes */
        for (int i = 0; i < nnodes; i++) {
            if (visited[i]) continue;

            /* iterate over saved edges */
            for (int j = 0; j < visnodes; j++) {
                double extra_milage =
                    dist(i, child->edges[j].i, inst) +
                    dist(i, child->edges[j].j, inst) -
                    dist(child->edges[j].i, child->edges[j].j, inst);

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
        int oldi = child->edges[edgei].i;
        int oldj = child->edges[edgei].j;
        child->edges[edgei] = (edge){oldi, next};
        child->edges[visnodes] = (edge){oldj, next};
        visnodes++;

        visited[next] = 1;

        child->zstar += best_extra_milage;

        if (VERBOSE) {
            printf("[VERBOSE] next: %d, break (%d, %d)\n", next + 1, oldi + 1,
                   oldj + 1);
        }
    }

    /* refine it! */
    if ((succ = edges_tosucc(child->edges, nnodes)) == NULL) {
        print_error("solver didnt produced a tour");
    }

    if (nnodes < 750)
        child->zstar += twoopt_refinement_notimelim(inst, succ, nnodes);

    /* store the succ as usual edges array */
    for (int i = 0; i < nnodes; i++) {
        child->edges[i] = (edge){i, succ[i]};
    }
    free(succ);
    free(chromosome);
    free(visited);

    return child;
}

solution TSPrandom(instance inst) {
    int nnodes = inst->nnodes;

    solution sol = create_solution(inst, GENETIC, nnodes);

    int* tour = (int*)calloc(nnodes, sizeof(int));
    int* succ = (int*)calloc(nnodes, sizeof(int));

    int permutation_steps = 2 * nnodes;

    for (int i = 0; i < nnodes; i++) tour[i] = i;

    for (int i = 0; i < permutation_steps; i++) {
        swap(&tour[rand() % nnodes], &tour[rand() % nnodes]);
    }

    for (int i = 0; i < nnodes; i++) {
        sol->zstar += dist(tour[i], tour[(i + 1) % nnodes], inst);
        succ[tour[i]] = tour[(i + 1) % nnodes];
    }

    if (nnodes < 200) {
        sol->zstar += twoopt_refinement_notimelim(inst, succ, nnodes);
    }
    /* store the succ as usual edges array */
    for (int j = 0; j < nnodes; j++) {
        sol->edges[j] = (edge){j, succ[j]};
    }
    free(succ);
    free(tour);

    return sol;
}
