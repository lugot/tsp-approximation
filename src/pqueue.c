#include "../include/pqueue.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../include/globals.h"

/* kinship helpers */
int parent(int i);
int left(int i);
int right(int i);

/* node helpers */
pqueue_node create_node(double key, int val);
void swap_node(pqueue_node* n, pqueue_node* m);

/* heap shift operations */
void shiftup(pqueue pq, int i);
void shiftdown(pqueue pq, int i);

/* kinship helpers */
int parent(int i) { return (i - 1) / 2; }
int left(int i) { return 2 * i + 1; }
int right(int i) { return 2 * i + 2; }

/* node helpers */
pqueue_node create_node(double key, int val) {
    pqueue_node n = (pqueue_node)calloc(1, sizeof(struct pqueue_node_t));
    n->key = key;
    n->val = val;

    return n;
}

void swap_node(pqueue_node* n, pqueue_node* m) {
    pqueue_node temp = *n;
    *n = *m;
    *m = temp;
}

/* heap shift operations */
void shiftup(pqueue pq, int i) {
    switch (pq->mode) {
        case MIN_HEAP: {
            while (i > 0 &&
                   pq->data[parent(i)]->key - pq->data[i]->key > -EPSILON) {
                swap_node(&pq->data[parent(i)], &pq->data[i]);
                i = parent(i);
            }
            break;
        }
        case MAX_HEAP: {
            while (i > 0 &&
                   pq->data[parent(i)]->key - pq->data[i]->key < EPSILON) {
                swap_node(&pq->data[parent(i)], &pq->data[i]);
                i = parent(i);
            }
            break;
        }
    }
}

void shiftdown(pqueue pq, int i) {
    switch (pq->mode) {
        case MIN_HEAP: {
            int min_index = i;

            int l = left(i);
            if (l <= pq->size &&
                pq->data[l]->key - pq->data[min_index]->key < EPSILON)
                min_index = l;

            int r = right(i);
            if (r <= pq->size &&
                pq->data[r]->key - pq->data[min_index]->key < EPSILON)
                min_index = r;

            if (i != min_index) {
                swap_node(&pq->data[i], &pq->data[min_index]);
                shiftdown(pq, min_index);
            }
            break;
        }
        case MAX_HEAP: {
            int max_index = i;

            int l = left(i);
            if (l <= pq->size &&
                pq->data[l]->key - pq->data[max_index]->key > -EPSILON)
                max_index = l;

            int r = right(i);
            if (r <= pq->size &&
                pq->data[r]->key - pq->data[max_index]->key > -EPSILON)
                max_index = r;

            if (i != max_index) {
                swap_node(&pq->data[i], &pq->data[max_index]);
                shiftdown(pq, max_index);
            }
            break;
        }
    }
}

/* priority queue operations */
pqueue pqueue_create(enum modes mode) {
    pqueue pq = (pqueue)calloc(1, sizeof(struct pqueue_t));

    pq->size = 0;
    pq->capacity = 10;
    pq->data = (pqueue_node*)calloc(10, sizeof(struct pqueue_node_t));
    pq->mode = mode;

    return pq;
}

int pqueue_empty(pqueue pq) { return pq->size == 0; }

int pqueue_top(pqueue pq) { return pq->data[0]->val; }

double pqueue_top_key(pqueue pq) { return pq->data[0]->key; }

int pqueue_pop(pqueue pq) {
    assert(pq->size != 0);

    int result = pq->data[0]->val;
    free(pq->data[0]);
    pq->data[0] = pq->data[pq->size - 1];

    pq->size--;

    shiftdown(pq, 0);

    return result;
}

void pqueue_push(pqueue pq, double key, int val) {
    pq->size++;

    if (pq->size >= pq->capacity) {
        pq->capacity = 2 * pq->capacity;
        pq->data = (pqueue_node*)realloc(
            pq->data, pq->capacity * sizeof(struct pqueue_node_t));
    }
    pq->data[pq->size - 1] = create_node(key, val);

    shiftup(pq, pq->size - 1);
}

void pqueue_print(pqueue pq) {
    if (pqueue_empty(pq)) {
        printf("<empty queue>\n");
        return;
    }

    pqueue temp = pqueue_create(pq->mode);

    while (!pqueue_empty(pq)) {
        double key = pq->data[0]->key;
        int val = pqueue_pop(pq);
        printf(" (%.3lf: , %d)", key, val);

        pqueue_push(temp, key, val);
    }
    while (!pqueue_empty(temp)) {
        double key = temp->data[0]->key;
        int val = pqueue_pop(temp);
        pqueue_push(pq, key, val);
    }
    printf("\n");
}
void pqueue_free(pqueue pq) {
    for (int i = 0; i < pq->size; i++) free(pq->data[i]);
    free(pq->data);
    free(pq);
}

/* topkqueue wrapper */
topkqueue topkqueue_create(int k) {
    topkqueue tk = (topkqueue)malloc(sizeof(struct topkqueue_t));
    tk->pq = pqueue_create(MAX_HEAP);
    tk->k = k;

    return tk;
}

void topkqueue_push(topkqueue tk, double key, int val) {
    if (tk->pq->size < tk->k) {
        /* if not enought pairs simply insert */
        pqueue_push(tk->pq, key, val);
    } else {
        /* new pair can enter if its key is lower than bigger one */
        if (key < pqueue_top_key(tk->pq)) {
            /* remove last pair (the one with bigger key) */
            pqueue_pop(tk->pq);
            /* insert newone */
            pqueue_push(tk->pq, key, val);
        }
    }
}

int topkqueue_randompick(topkqueue tk) {
    unsigned int seedp = time(NULL);
    int npops = rand_r(&seedp) % tk->pq->size;
    while (npops--) pqueue_pop(tk->pq);
    int ans = pqueue_top(tk->pq);

    /* let's also free the nodes, there the overhead of shiftdown only */
    while (!pqueue_empty(tk->pq)) pqueue_pop(tk->pq);

    return ans;
}

void topkqueue_print(topkqueue tk) { pqueue_print(tk->pq); }

void topkqueue_free(topkqueue tk) {
    pqueue_free(tk->pq);
    free(tk);
}
