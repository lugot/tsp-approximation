#include "../include/adjlist.h"

#include <assert.h>
#include <stdlib.h>

adjlist adjlist_create(int N) {
    adjlist l = (adjlist)calloc(1, sizeof(struct adjlist_t));

    l->pairs = (pair*)calloc(N, sizeof(struct pair_t));
    for (int i = 0; i < N; i++) {
        l->pairs[i] = (pair)calloc(1, sizeof(struct pair_t));
        l->pairs[i]->a = l->pairs[i]->b = -1;
    }

    l->npairs = N;

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

pair* adjlist_loose_ends(adjlist l, int* npairs) {
    int N = l->npairs;
    pair* pairs = l->pairs;

    pair* ans = (pair*)calloc(N, sizeof(struct pair_t));
    *npairs = 0;

    int front, back;
    int* q = (int*)malloc(N * sizeof(int));
    int* visited = (int*)calloc(N, sizeof(int));
    for (int i = 0; i < N; i++) q[i] = -1;

    for (int i = 0; i < N; i++) {
        if (visited[i]) continue;

        front = back = 0;
        q[++back] = i;

        pair p = (pair)calloc(1, sizeof(struct pair_t));
        p->a = p->b = -1;

        while (front < back) {
            int act = q[front++];

            /* add to the solution */
            if (pairs[act]->a == -1 || pairs[act]->b == -1) {
                if (p->a == -1) {
                    p->a = act;
                } else {
                    p->b = act;
                }
            }

            /* add to the queue */
            if (pairs[act]->a != -1 && !visited[act]) {
                visited[act] = 1;
                q[++back] = pairs[act]->b;
            }
            if (pairs[act]->b != -1 && !visited[act]) {
                visited[act] = 1;
                q[++back] = pairs[act]->b;
            }
        }

        /* add to the solution */
        ans[*npairs] = p;
        (*npairs)++;
        free(p);
    }
    free(q);
    free(visited);

    return ans;
}
