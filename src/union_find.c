#include "../include/union_find.h"

#include <stdlib.h>
#include <string.h>

#include "../include/utils.h"

union_find uf_create(int N) {
    union_find uf = (union_find)calloc(1, sizeof(struct union_find_t));

    uf->p = (int*)calloc(N, sizeof(int));
    uf->rank = (int*)calloc(N, sizeof(int));
    uf->size_of_set = (int*)calloc(N, sizeof(int));
    uf->next = (int*)calloc(N, sizeof(int));

    /* initialize data structure */
    for (int i = 0; i < N; i++) uf->p[i] = i;
    /* rank already initialized by calloc */
    /*memset(uf->size_of_set, 1, N*sizeof(int));*/
    for (int i = 0; i < N; i++) uf->size_of_set[i] = 1;
    for (int i = 0; i < N; i++) uf->next[i] = i;

    uf->num_sets = N;

    return uf;
}

int uf_find_set(union_find uf, int i) {
    return (uf->p[i] == i) ? i : (uf->p[i] = uf_find_set(uf, uf->p[i]));
}

int uf_same_set(union_find uf, int i, int j) {
    return uf_find_set(uf, i) == uf_find_set(uf, j);
}

int uf_set_size(union_find uf, int i) {
    return uf->size_of_set[uf_find_set(uf, i)];
}

void uf_union_set(union_find uf, int i, int j) {
    if (!uf_same_set(uf, i, j)) {
        int x = uf_find_set(uf, i);
        int y = uf_find_set(uf, j);

        if (uf->rank[x] > uf->rank[y]) swap(&x, &y);
        uf->p[x] = y;
        if (uf->rank[x] == uf->rank[y]) uf->rank[y]++;
        uf->size_of_set[y] += uf->size_of_set[x];
        uf->num_sets--;

        swap(&uf->next[x], &uf->next[y]);
    }
}

void uf_free(union_find uf) {
    free(uf->p);
    free(uf->rank);
    free(uf->size_of_set);

    free(uf);
}
