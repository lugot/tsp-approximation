#ifndef _PQUEUE_H_
#define _PQUEUE_H_

#include "event.h"
#include <stddef.h>

typedef struct pqueue_node_t {
    double key;
    event ev;
} *pqueue_node; 

typedef struct pqueue_t {
	struct pqueue_node_t** data;
	size_t size;
	size_t capacity;
} *pqueue;


/* priority queue operations */
pqueue pqueue_create();
int pqueue_empty(pqueue pq);
event pqueue_top(pqueue pq);
event pqueue_pop(pqueue pq);
void pqueue_push(pqueue pq, double key, event ev);
void pqueue_print(pqueue pq);

#endif /* _PQUEUE_H_ */
