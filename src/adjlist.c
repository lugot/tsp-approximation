#include "../include/adjlist.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

adjlist adjlist_create(int N) {
    adjlist l = (adjlist)calloc(1, sizeof(struct adjlist_t));

    l->neigh = (int(*)[2])malloc(N * sizeof(int[2]));
    l->N = N;

    l->visited = (int*)calloc(N, sizeof(int));

    adjlist_hard_reset(l);

    return l;
}
void adjlist_free(adjlist l) {
    free(l->neigh);
    free(l->visited);

    free(l);
}

void adjlist_add_edge(adjlist l, int i, int j) {
    int(*neigh)[2] = l->neigh;

    assert(neigh[i][0] == -1 || neigh[i][1] == -1);
    assert(neigh[j][0] == -1 || neigh[j][1] == -1);

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
    int(*neigh)[2] = l->neigh;
    int N = l->N;

    for (int i = 0; i < N; i++) {
        printf("%d: pair %d, %d\n", i + 1, neigh[i][0] + 1, neigh[i][1] + 1);
    }
}
void adjlist_reset(adjlist l) {
    /* reset state: visited vector, i and bit */
    memset(l->visited, 0, l->N * sizeof(int));
    l->i = 0;
    l->bit = 0;
}
void adjlist_hard_reset(adjlist l) {
    /* reset state and inserted arcs */
    adjlist_reset(l);

    memset(l->neigh, -1, l->N * sizeof(int[2]));
}

int adjlist_get_edge(adjlist l, int* u, int* v) {
    /* writing l-> is so boring.. */
    int(*neigh)[2] = l->neigh;
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

int adjlist_get_loose_ends(adjlist l, int* end1, int* end2) {
    /* writing l-> is so boring.. */
    int(*neigh)[2] = l->neigh;
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

int* adjlist_get_subtour(adjlist l, int* subsize) {
    /* writing l-> is so boring.. */
    int(*neigh)[2] = l->neigh;
    int N = l->N;
    int* visited = l->visited;

    int front, back;
    int* q = (int*)malloc(N * sizeof(int));
    q = (int*)memset(q, -1, N * sizeof(int));

    int* subtour = (int*)malloc(N * sizeof(int));

    while (l->i < N && (visited[l->i] || (neigh[l->i][0] == neigh[l->i][1])))
        l->i++;
    if (l->i == N) {
        *subsize = 0;
        free(subtour);
        free(q);
        return NULL;
    }

    front = back = 0;
    q[back++] = l->i;
    visited[l->i] = 1;
    *subsize = 0;
    subtour[(*subsize)++] = l->i;

    while (front < back) {
        int act = q[front++];
        int u, v;
        u = neigh[act][0];
        v = neigh[act][1];

        /* add to the queue and to the solution*/
        if (u != -1 && !visited[u]) {
            visited[u] = 1;
            q[back++] = u;
            subtour[(*subsize)++] = u;
        }
        if (v != -1 && !visited[v]) {
            visited[v] = 1;
            q[back++] = v;
            subtour[(*subsize)++] = v;
        }
    }

    free(q);

    if (*subsize == N) {
        free(subtour);
        return NULL;
    } else {
        return subtour;
    }
}

int adjlist_single_tour(adjlist l) {
    int subsize;
    int* subtour = adjlist_get_subtour(l, &subsize);
    free(subtour);

    adjlist_reset(l);

    return subsize == l->N;
}
