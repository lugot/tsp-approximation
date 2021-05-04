#include "../include/adjlist.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

adjlist adjlist_create(int N) {
    adjlist l = (adjlist)calloc(1, sizeof(struct adjlist_t));

    l->neigh = (int**)calloc(N, sizeof(int*));
    for (int i = 0; i < N; i++) {
        l->neigh[i] = (int*)calloc(2, sizeof(int));
        l->neigh[i][0] = l->neigh[i][1] = -1;
    }

    l->N = N;

    l->visited = (int*)calloc(N, sizeof(int));
    l->i = 0;
    l->bit = 0;

    return l;
}
void adjlist_free(adjlist l) {
    for (int i = 0; i < l->N; i++) free(l->neigh[i]);
    free(l->visited);

    free(l);
}

void adjlist_add_arc(adjlist l, int i, int j) {
    int** neigh = l->neigh;

    assert(neigh[i][0] == -1 || neigh[i][1] == -1);
    assert(neigh[j][1] == -1 || neigh[j][1] == -1);

    if (neigh[i][0] == -1) {
        neigh[i][0] = j;
    } else {
        neigh[i][1] = j;
    }

    if (neigh[j][0] == -1) {
        neigh[j][0] = i;
    } else {
        neigh[j][1] = i;
    }
}
void adjlist_print(adjlist l) {
    int** neigh = l->neigh;
    int N = l->N;

    for (int i = 0; i < N; i++) {
        printf("%d: pair %d, %d\n", i + 1, neigh[i][0] + 1, neigh[i][1] + 1);
    }
}
void adjlist_reset(adjlist l) {
    l->visited = (int*)memset(l->visited, 0, l->N * sizeof(int));
    l->i = 0;
    l->bit = 0;
}

int adjlist_get_arc(adjlist l, int* u, int* v) {
    /* writing l-> is so boring.. */
    int** neigh = l->neigh;
    int N = l->N;
    int arc_found = 0;

    do {
        while (l->i < N && (neigh[l->i][0] == -1)) l->i++;
        if (l->i == N) return 0;

        if (neigh[l->i][l->bit] != -1) {
            *u = l->i;
            *v = neigh[l->i][l->bit];

            if (*u < *v) arc_found = 1;
        }
        l->bit = (l->bit + 1) % 2;
        if (l->bit == 0) l->i++;

    } while (!arc_found);

    return 1;
}

int adjlist_loose_ends(adjlist l, int* end1, int* end2) {
    /* writing l-> is so boring.. */
    int** neigh = l->neigh;
    int N = l->N;
    int* visited = l->visited;

    int front, back;
    int* q = (int*)malloc(N * sizeof(int));
    q = (int*)memset(q, -1, N * sizeof(int));

    while (l->i < N && (visited[l->i] || (neigh[l->i][0] == neigh[l->i][1])))
        l->i++;
    if (l->i == N) return 0;

    front = back = 0;
    q[back++] = l->i;
    visited[l->i] = 1;

    *end1 = *end2 = -1;
    while (front < back) {
        int act = q[front++];
        int u, v;
        u = neigh[act][0];
        v = neigh[act][1];

        /* add to the solution */
        if (u == -1 || v == -1) {
            if (*end1 == -1) {
                *end1 = act;
            } else {
                *end2 = act;
            }
        }

        /* add to the queue */
        if (u != -1 && !visited[u]) {
            visited[u] = 1;
            q[back++] = u;
        }
        if (v != -1 && !visited[v]) {
            visited[v] = 1;
            q[back++] = v;
        }
    }
    free(q);

    return 1;
}

/* pair* adjlist_loose_ends(adjlist l, int* nneigh) { */
/*     int N = l->nneigh; */
/*     pair* neigh = l->pairs; */

/*     pair* ans = (pair*)calloc(N, sizeof(struct pair_t)); */
/*     *nneigh = 0; */

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
/*         p[0] = p[1] = -1; */

/*         int skip = 0; */

/*         while (front < back) { */
/*             int act = q[front++]; */
/*             int i, j; */
/*             i = neigh[act][0]; */
/*             j = neigh[act][1]; */
/*             printf("act: %d, i: %d, j: %d\n", act + 1, i + 1, j + 1); */

/*             if (i == -1 && j == -1) { */
/*                 skip = 1; */
/*                 break; */
/*             } */

/*             /1* add to the solution *1/ */
/*             if (i == -1 || j == -1) { */
/*                 if (p[0] == -1) { */
/*                     p[0] = act; */
/*                 } else { */
/*                     p[1] = act; */
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
/*         if (!skip) ans[(*nneigh)++] = p; */
/*     } */
/*     free(q); */
/*     free(visited); */

/*     return ans; */
/* } */
