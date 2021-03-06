#include "../include/metaheuristics.h"

#include <assert.h>
#include <float.h>
#include <time.h>
#include <unistd.h>

#include "../include/constructives.h"
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

            tracker_add(sol->t, stopwatch(&s, &e), sol->zstar);

            /* restart k from initial value */
            k = VNS_K_START;

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
    solution sol = create_solution(inst, TABU_SEACH_RANDOM, nnodes);
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

    srand(time(NULL));

    /* initialize total wall-clock time of execution time */
    struct timespec s, e;
    s.tv_sec = e.tv_sec = -1;
    stopwatch(&s, &e);

    int nnodes = inst->nnodes;
    solution sol = create_solution(inst, GENETIC, nnodes);
    sol->zstar = DBL_MAX;
    solution champion;

    solution* population = (solution*)calloc(GENETIC_N, sizeof(solution));
    solution* offspring = (solution*)calloc(GENETIC_K, sizeof(solution));

    /* generate initial populaiton */
    for (int i = 0; i < GENETIC_N; i++) {
        if ((double)rand() / (double)RAND_MAX < GENETIC_PERC_LUCKYCHILDREN) {
            population[i] = TSPgrasp(inst, 1);
            population[i]->model_type = GENETIC;
        } else {
            population[i] = TSPrandom(inst);
        }
    }
    if (VERBOSE) printf("[VERBOSE] done generating population\n");

    /* update timelimit */
    double timelimit = inst->params->timelimit;
    inst->params->timelimit = timelimit *= GENETIC_PERC_TIME;
    while (stopwatch(&s, &e) / 1000.0 < inst->params->timelimit) {
        /* pick randomly two parents */
        for (int i = 0; i < GENETIC_K; i++) {
            int parent1_idx, parent2_idx;

            parent1_idx = rand() % GENETIC_N;
            parent2_idx = rand() % GENETIC_N;
            while (parent2_idx == parent1_idx) parent2_idx = rand() % GENETIC_N;

            offspring[i] =
                child(inst, population[parent1_idx], population[parent2_idx]);
        }
        if (stopwatch(&s, &e) / 1000.0 > inst->params->timelimit) {
            break;
        }
        if (VERBOSE) printf("[VERBOSE] done generating children\n");

        /* insert randomly the offspring in the population and kill parents but
         * keep the best ones*/
        int* best_sols = (int*)calloc(GENETIC_N, sizeof(int));
        pair* objs = (pair*)calloc(GENETIC_N, sizeof(pair));
        for (int i = 0; i < GENETIC_N; i++) {
            objs[i] = (pair){population[i]->zstar, i};
        }
        qsort(objs, GENETIC_N, sizeof(pair), paircmp);

        for (int i = 0; i < GENETIC_NBESTSOLS; i++) {
            best_sols[objs[i].x] = 1;
        }

        for (int i = 0; i < GENETIC_K; i++) {
            int parent = rand() % GENETIC_N;
            while (best_sols[parent]) parent = rand() % GENETIC_N;

            free_solution(population[parent]);
            population[parent] = offspring[i];
        }
        if (VERBOSE) printf("[VERBOSE] done killing\n");

        /* mutate the population */
        int timelimit_reached = 0;
        for (int i = 0; i < GENETIC_NMUTATIONS; i++) {
            if (timelimit_reached) break;

            int tomutate = rand() % GENETIC_N;
            int* succ;
            if ((succ = edges_tosucc(population[tomutate]->edges, nnodes)) ==
                NULL) {
                print_error("no solution found!");
            }

            for (int j = 0; j < GENETIC_NMOVES; j++) {
                int a, b;
                double delta;
                a = b = 0;

                delta = twoopt_pick(inst, succ, &a, &b);
                if (stopwatch(&s, &e) / 1000.0 > inst->params->timelimit) {
                    timelimit_reached = 1;
                    break;
                }
                if (delta == 0.0) break; /* no more moves possible */

                if (EXTRA_VERBOSE) {
                    printf("[VERBOSE] refinement on %d, %d delta %lf\n", a, b,
                           delta);
                }

                /* actually perform the move */
                twoopt_move(succ, nnodes, a, b);
                /* update the objective */
                population[tomutate]->zstar += delta;
            }

            for (int i = 0; i < nnodes; i++) {
                population[tomutate]->edges[i] = (edge){i, succ[i]};
            }

            free(succ);
        }
        if (VERBOSE) printf("[VERBOSE] done mutating\n");

        /* track champion */
        double best_zstar = DBL_MAX;
        for (int i = 0; i < GENETIC_N; i++) {
            if (population[i]->zstar < best_zstar) {
                best_zstar = population[i]->zstar;

                champion = population[i];
            }
        }

        if (VERBOSE) {
            printf("[VERBOSE] champion zstar: %lf\n", champion->zstar);
        }

        if (champion->zstar < sol->zstar) {
            sol->zstar = champion->zstar;
            tracker_add(sol->t, stopwatch(&s, &e), sol->zstar);
        }
    }

    /* last champion becames the solution returned */
    sol->zstar = champion->zstar;
    for (int i = 0; i < nnodes; i++) {
        sol->edges[i] = (edge){champion->edges[i].i, champion->edges[i].j};
    }

    /* refine best solution if we have time! */
    int* succ = (int*)calloc(nnodes, sizeof(int));
    if ((succ = edges_tosucc(sol->edges, nnodes)) == NULL) {
        print_error("solver didnt produced a tour");
    }

    /* reset times */
    inst->params->timelimit = timelimit * (1 - GENETIC_PERC_TIME);
    s.tv_sec = e.tv_sec = -1;
    stopwatch(&s, &e);

    sol->zstar += twoopt_refinement(inst, succ, nnodes, &s, &e);
    free(succ);

    for (int i = 0; i < GENETIC_N; i++) {
        free_solution(population[i]);
    }
    free(population);
    free(offspring);

    return sol;
}

solution child(instance inst, solution parent1, solution parent2) {
    if (VERBOSE) {
        printf("[VERBOSE] creating children!\n\n");
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

        if (EXTRA_VERBOSE) {
            printf("[VERBOSE] next: %d, break (%d, %d)\n", next + 1, oldi + 1,
                   oldj + 1);
        }
    }

    /* refine it! */
    if ((succ = edges_tosucc(child->edges, nnodes)) == NULL) {
        print_error("solver didnt produced a tour");
    }

    /* if (nnodes < 750) */
    /*     child->zstar += twoopt_refinement_notimelim(inst, succ, nnodes);
     */

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

    unsigned int seedp = rand();
    solution sol = create_solution(inst, GENETIC, nnodes);

    int* succ = randomtour(nnodes, seedp);
    for (int i = 0; i < nnodes; i++) {
        sol->edges[i] = (edge){i, succ[i]};
        sol->zstar += dist(i, succ[i], inst);
    }
    free(succ);

    return sol;
}

/* OLD TSPrandom */
/* int* tour = (int*)calloc(nnodes, sizeof(int)); */
/* int* succ = (int*)calloc(nnodes, sizeof(int)); */
/*  */
/* int permutation_steps = 2 * nnodes; */
/*  */
/* for (int i = 0; i < nnodes; i++) tour[i] = i; */
/*  */
/* for (int i = 0; i < permutation_steps; i++) { */
/*     swap(&tour[rand() % nnodes], &tour[rand() % nnodes]); */
/* } */
/*  */
/* for (int i = 0; i < nnodes; i++) { */
/*     sol->zstar += dist(tour[i], tour[(i + 1) % nnodes], inst); */
/*     succ[tour[i]] = tour[(i + 1) % nnodes]; */
/* } */
/*  */
/* #<{(| store the succ as usual edges array |)}># */
/* for (int j = 0; j < nnodes; j++) { */
/*     sol->edges[j] = (edge){j, succ[j]}; */
/* } */
/* free(succ); */
/* free(tour); */
