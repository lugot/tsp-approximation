#include "../include/union_find.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "../include/utils.h"

/* union find members */
union_find uf_create(int N) {
    union_find uf = (union_find)calloc(1, sizeof(struct union_find_t));

    /* initialize data structure */
    uf->p = (int*)calloc(N, sizeof(int));
    for (int i = 0; i < N; i++) uf->p[i] = i;
    uf->rank = (int*)calloc(N, sizeof(int));
    uf->size_of_set = (int*)calloc(N, sizeof(int));
    memset(uf->size_of_set, 1, N * sizeof(int));

    uf->N = uf->nsets = N;

    uf->tns = (tree_node*)malloc(N * sizeof(struct tree_node_t));
    for (int i = 0; i < N; i++) uf->tns[i] = tn_create();
    uf->stack = (int*)calloc(N, sizeof(int));
    uf->si = -1;
    uf->visited = (int*)calloc(N, sizeof(int));

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
        tn_add_son(uf->tns[i], j);
        tn_add_son(uf->tns[j], i);

        int x = uf_find_set(uf, i);
        int y = uf_find_set(uf, j);

        if (uf->rank[x] > uf->rank[y]) swap(&x, &y);

        uf->p[x] = y;

        if (uf->rank[x] == uf->rank[y]) uf->rank[y]++;
        uf->size_of_set[y] += uf->size_of_set[x];
        uf->nsets--;
    }
}
int uf_postorder(union_find uf, int* next) {
    assert(uf->nsets == 1);

    /* start the visit */
    if (uf->si == -1) {
        uf->stack[0] = uf_find_set(uf, 0);
        uf->si = 1;
        uf->visited[uf->stack[0]] = 1;
    }

    /* end visit */
    if (uf->si == 0) return 0;

    *next = uf->stack[uf->si - 1]; /* pop */
    uf->si--;
    uf->visited[*next] = 1;

    /* do not even ask */
    for (int i = 0; i < uf->tns[*next]->size; i++) {
        int son = uf->tns[*next]->sons[i];

        if (uf->visited[son]) continue;

        uf->stack[uf->si++] = son;
    }

    return 1;
}

void uf_free(union_find uf) {
    free(uf->p);
    free(uf->rank);
    free(uf->size_of_set);
    for (int i = 0; i < uf->N; i++) tn_free(uf->tns[i]);
    free(uf->tns);
    free(uf->stack);
    free(uf->visited);

    free(uf);
}

/* tree node members */
tree_node tn_create() {
    tree_node tn = (tree_node)malloc(sizeof(struct tree_node_t));

    tn->size = 0;
    tn->capacity = 10;
    tn->sons = (int*)calloc(tn->capacity, sizeof(int));

    return tn;
}
void tn_add_son(tree_node tn, int s) {
    tn->size++;

    if (tn->size >= tn->capacity) {
        tn->capacity = 2 * tn->capacity;
        tn->sons = (int*)realloc(tn->sons, tn->capacity * sizeof(int));
    }
    tn->sons[tn->size - 1] = s;
}
void tn_free(tree_node tn) {
    free(tn->sons);
    free(tn);
}
