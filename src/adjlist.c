#include "../include/adjlist.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

adjlist adjlist_create(int N) {
    adjlist l = (adjlist)calloc(1, sizeof(struct adjlist_t));

    l->pairs = (pair*)calloc(N, sizeof(struct pair_t));
    for (int i = 0; i < N; i++) {
        l->pairs[i] = (pair)calloc(1, sizeof(struct pair_t));
        l->pairs[i]->a = l->pairs[i]->b = -1;
    }

    l->npairs = N;

    l->visited = (int*)calloc(N, sizeof(int));
    l->i = 0;

    return l;
}
void adjlist_free(adjlist l) {
    for (int i = 0; i < l->npairs; i++) free(l->pairs[i]);
    free(l);
}

void adjlist_add_arc(adjlist l, int i, int j) {
    pair* pairs = l->pairs;

    assert(pairs[i]->a == -1 || pairs[i]->b == -1);
    assert(pairs[j]->a == -1 || pairs[j]->b == -1);

    if (pairs[i]->a == -1) {
        pairs[i]->a = j;
    } else {
        pairs[i]->b = j;
    }

    if (pairs[j]->a == -1) {
        pairs[j]->a = i;
    } else {
        pairs[j]->b = i;
    }
}
void adjlist_print(adjlist l) {
    int N = l->npairs;
    pair* pairs = l->pairs;

    for (int i = 0; i < N; i++) {
        printf("pair %d, %d\n", pairs[i]->a + 1, pairs[i]->b + 1);
    }
}
void adjlist_reset(adjlist l) {
    l->visited = (int*)memset(l->visited, 0, l->npairs * sizeof(int));
    l->i = 0;
}

int adjlist_loose_ends(adjlist l, int* end1, int* end2) {
    int N = l->npairs;
    pair* pairs = l->pairs;

    int front, back;
    int* q = (int*)malloc(N * sizeof(int));
    q = (int*)memset(q, -1, N * sizeof(int));

    while (l->i < N &&
           (l->visited[l->i] || (l->pairs[l->i]->a == l->pairs[l->i]->b)))
        l->i++;
    if (l->i == N) return 0;

    front = back = 0;
    q[back++] = l->i;
    l->visited[l->i] = 1;

    *end1 = *end2 = -1;
    while (front < back) {
        int act = q[front++];
        int i, j;
        i = pairs[act]->a;
        j = pairs[act]->b;

        /* add to the solution */
        if (i == -1 || j == -1) {
            if (*end1 == -1) {
                *end1 = act;
            } else {
                *end2 = act;
            }
        }

        /* add to the queue */
        if (i != -1 && !l->visited[i]) {
            l->visited[i] = 1;
            q[back++] = i;
        }
        if (j != -1 && !l->visited[j]) {
            l->visited[j] = 1;
            q[back++] = j;
        }
    }
    free(q);

    return 1;
}

/* pair* adjlist_loose_ends(adjlist l, int* npairs) { */
/*     int N = l->npairs; */
/*     pair* pairs = l->pairs; */

/*     pair* ans = (pair*)calloc(N, sizeof(struct pair_t)); */
/*     *npairs = 0; */

/*     int front, back; */
/*     int* q = (int*)malloc(N * sizeof(int)); */
/*     int* visited = (int*)calloc(N, sizeof(int)); */
/*     for (int i = 0; i < N; i++) q[i] = -1; */

/*     for (int k = 0; k < N; k++) { */
/*         if (visited[k]) continue; */

/*         front = back = 0; */
/*         q[back++] = k; */
/*         visited[k] = 1; */

/*         pair p = (pair)calloc(1, sizeof(struct pair_t)); */
/*         p->a = p->b = -1; */

/*         int skip = 0; */

/*         while (front < back) { */
/*             int act = q[front++]; */
/*             int i, j; */
/*             i = pairs[act]->a; */
/*             j = pairs[act]->b; */
/*             printf("act: %d, i: %d, j: %d\n", act + 1, i + 1, j + 1); */

/*             if (i == -1 && j == -1) { */
/*                 skip = 1; */
/*                 break; */
/*             } */

/*             /1* add to the solution *1/ */
/*             if (i == -1 || j == -1) { */
/*                 if (p->a == -1) { */
/*                     p->a = act; */
/*                 } else { */
/*                     p->b = act; */
/*                 } */
/*             } */

/*             /1* add to the queue *1/ */
/*             if (i != -1 && !visited[i]) { */
/*                 visited[i] = 1; */
/*                 q[back++] = i; */
/*             } */
/*             if (j != -1 && !visited[j]) { */
/*                 visited[j] = 1; */
/*                 q[back++] = j; */
/*             } */
/*         } */

/*         /1* add to the solution *1/ */
/*         if (!skip) ans[(*npairs)++] = p; */
/*     } */
/*     free(q); */
/*     free(visited); */

/*     return ans; */
/* } */
