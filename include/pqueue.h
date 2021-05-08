#ifndef INCLUDE_PQUEUE_H_
#define INCLUDE_PQUEUE_H_

#include <stddef.h>

enum modes { MAX_HEAP, MIN_HEAP };

typedef struct pqueue_node_t {
    double key;
    int val;
} * pqueue_node;

typedef struct pqueue_t {
    struct pqueue_node_t** data;
    int size;
    int capacity;
    enum modes mode;
} * pqueue;

typedef struct topkqueue_t {
    pqueue pq;
    int k;
} * topkqueue;

/* priority queue operations */
pqueue pqueue_create(enum modes mode);
int pqueue_empty(pqueue pq);
int pqueue_top(pqueue pq);
double pqueue_top_key(pqueue pq);
int pqueue_pop(pqueue pq);
void pqueue_push(pqueue pq, double key, int val);
void pqueue_print(pqueue pq);
void pqueue_free(pqueue pq);

/* topkqueue: wrapper of priority list which mantains only the top-k elements */
topkqueue topkqueue_create(int k);
void topkqueue_push(topkqueue tk, double key, int val);
int topkqueue_randompick(topkqueue tk);
void topkqueue_print(topkqueue tk);
void topkqueue_free(topkqueue tk);


#endif  // INCLUDE_PQUEUE_H_
