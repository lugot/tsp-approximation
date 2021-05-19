#define _GNU_SOURCE

#include <assert.h>
#include <stdlib.h>
#include <time.h>

#include "../include/globals.h"
#include "../include/refinements.h"
#include "../include/utils.h"
#include "../include/union_find.h"

double twoopt_pick(instance inst, int* succ, int* a, int* b);
void twoopt_move(int* succ, int nnodes, int a, int b);
double threeopt_refinement(instance inst, int* succ, int nnodes);
double threeopt_pick(instance inst, int* succ, int* a, int* b, int* c);
void threeopt_move(int* succ, int nnodes, int a, int b, int c, instance inst);

double twoopt_refinement(instance inst, int* succ, int nnodes) {
    double improvement = 0.0;

    int a, b;
    double delta = 1.0; /* using 0.0 as sentinel */
    a = b = 0;

    /* iterate over 2opt moves until not improvable */
    while ((delta = twoopt_pick(inst, succ, &a, &b)) != 0.0) {
        if (EXTRA) {
            printf("[VERBOSE] refinement on %d, %d delta %lf\n", a, b, delta);
        }

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
    {}
}

double threeopt_refinement(instance inst, int* succ, int nnodes) {
    double improvement = 0.0;

    int a, b, c;
    double delta = 1.0; /* using 0.0 as sentinel */
    a = b = c = 0;

    /* iterate over 2opt moves until not improvable */
    while ((delta = threeopt_pick(inst, succ, &a, &b, &c)) != 0.0) {
        if (EXTRA) {
            printf("[VERBOSE] refinement on %d, %d, %d delta %lf\n", a, b, c,
                   delta);
        }

        /* actually perform the move */
        threeopt_move(succ, nnodes, a, b, c, inst);

        /* update the objective */
        improvement += delta; /* delta should be negative */
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
            for (int k = j + 1; k < nnodes; k++) {
                if (i == succ[j] || i == succ[k]) continue;
                if (j == succ[i] || j == succ[k]) continue;
                if (k == succ[i] || k == succ[j]) continue;

                double delta[4];
                for (int n = 0; n < 4; n++)
                    delta[n] =
                        -(dist(i, succ[i], inst) + dist(j, succ[j], inst) +
                          dist(k, succ[k], inst));

                int g = i;
                int swap = 0;
                int tj, tk;

                // makes k the third visited node
                while (succ[g] != i) {
                    if (swap == 0) {
                        if (succ[g] == k) swap = 1;
                        if (succ[g] == j) swap = -1;
                    }
                    g = succ[g];
                }

                if (swap == 1) {
                    tk = j;
                    tj = k;
                } else {
                    tj = j;
                    tk = k;
                }

                // 5 cases: preserving cases reducing to 2-opt and subtours
                delta[0] += dist(i, tj, inst) + dist(succ[i], tk, inst) +
                            dist(succ[tj], succ[tk],
                                 inst);  // i, j, succ[i], k, succ[j], succ[k]
                delta[1] += dist(i, succ[tj], inst) + dist(tk, succ[i], inst) +
                            dist(tj, succ[tk],
                                 inst);  // i, succ[j], k, succ[i], j, succ[k]
                delta[2] += dist(i, succ[tj], inst) + dist(tk, tj, inst) +
                            dist(succ[i], succ[tk],
                                 inst);  // i, succ[j], k, j, succ[i], succ[k]
                delta[3] += dist(i, tk, inst) + dist(succ[tj], succ[i], inst) +
                            dist(tj, succ[tk],
                                 inst);  // i, k, succ[j], succ[i], j, succ[k]

                for (int n = 0; n < 4; n++)
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

    delta[0] =
        dist(a, b, inst) + dist(succ[a], c, inst) +
        dist(succ[b], succ[c], inst);  // i, j, succ[i], k, succ[j], succ[k]
    delta[1] = dist(a, succ[b], inst) + dist(c, succ[a], inst) +
               dist(b, succ[c], inst);  // i, succ[j], k, succ[i], j, succ[k]
    delta[2] =
        dist(a, succ[b], inst) + dist(c, b, inst) +
        dist(succ[a], succ[c], inst);  // i, succ[j], k, j, succ[i], succ[k]
    delta[3] = dist(a, c, inst) + dist(succ[b], succ[a], inst) +
               dist(b, succ[c], inst);  // i, k, succ[j], succ[i], j, succ[k]

    double best = delta[0];

    for (int n = 1; n < 4; n++)
        if (delta[n] < best) {
            best = delta[n];
            index = n;
        }

    best = best - (dist(a, succ[a], inst) + dist(b, succ[b], inst) +
                   dist(c, succ[c], inst));

    // if necessary swap tour from aprime to b
    if (index == 0 || index == 2) {
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

    // if necessary swap tour from bprime to c
    if (index == 0 || index == 3) {
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

    if (index == 0) {
        succ[a] = b;
        succ[aprime] = c;
        succ[bprime] = cprime;
    }
    if (index == 1) {
        succ[a] = bprime;
        succ[c] = aprime;
        succ[b] = cprime;
    }
    if (index == 2) {
        succ[a] = bprime;
        succ[c] = b;
        succ[aprime] = cprime;
    }
    if (index == 3) {
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
    unsigned int seedp = time(NULL);

    /* phase 1: embedded successor array preparation */
    int embsize = 2 * strength; /* embedded size */

    int nselected = 0;
    int* selected = (int*)malloc(strength * sizeof(int));
    while (nselected < strength) {
        int candidate = rand_r(&seedp) % nnodes;

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

    /* *looking at the heap as being vito corleone* look how they massacred my
     * boy*/
    uf_free(uf);
    free(pathlenghts);
    free(selected);
    free(buffer);
    free(embsucc);
    free(map);
}
