#include "pqueue.h"
#include "globals.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

typedef pqueue_node node;

/* kinship helpers */
int parent(size_t i);
int left(size_t i);
int right(size_t i); 

/* node helpers */
node create_node(double key, event ev);
void swap_node(node* n, node* m);

/* heap shift operations */
void shiftup(pqueue pq, size_t i);
void shiftdown(pqueue pq, size_t i);


/* kinship helpers */
int parent(size_t i) {
    return (i-1)/2; 
}
int left(size_t i) { 
    return 2*i + 1; 
}
int right(size_t i) { 
    return 2*i + 2; 
}

/* node helpers */
node create_node(double key, event ev) {
    node n = (node) calloc(1, sizeof(struct pqueue_node_t));
    n->key = key;
    n->ev = ev;

    return n;
}
void swap_node(node* n, node* m) {
    node temp = *n;
	*n = *m;
	*m = temp;
}

/* heap shift operations */
void shiftup(pqueue pq, size_t i) {
    while (i > 0 && pq->data[parent(i)]->key - pq->data[i]->key < GLOBALS.EPS) {
        swap_node(&pq->data[parent(i)], &pq->data[i]); 
        i = parent(i); 
    } 
} 
void shiftdown(pqueue pq, size_t i) {
	size_t maxIndex = i; 
  
	size_t l = left(i);
	if (l <= pq->size && pq->data[l]->key - pq->data[maxIndex]->key > -GLOBALS.EPS)
		maxIndex = l;

	size_t r = right(i);
	if (r <= pq->size && pq->data[r]->key - pq->data[maxIndex]->key > -GLOBALS.EPS)
		maxIndex = r;

	if (i != maxIndex){
		swap_node(&pq->data[i], &pq->data[maxIndex]);
		shiftdown(pq, maxIndex);
	}
} 

/* priority queue operations */
pqueue pqueue_create() {
    pqueue pq = (pqueue) calloc(1, sizeof(struct pqueue_t));

	pq->size = 0;
	pq->capacity = 10;
	pq->data = (node*) calloc(10, sizeof(struct pqueue_node_t));

    return pq;
}
int pqueue_empty(pqueue pq) {
	return pq->size == 0;
}
event pqueue_top(pqueue pq) {
	return pq->data[0]->ev;
}
double pqueue_top_key(pqueue pq) {
    return pq->data[0]->key;
}
event pqueue_pop(pqueue pq) {
    assert(pq->size != 0);

	event result = pq->data[0]->ev; 
	pq->data[0] = pq->data[pq->size-1];
	pq->size--;
  
	shiftdown(pq, 0);

	return result; 
} 
void pqueue_push(pqueue pq, double key, event ev) {
	pq->size++;

	if (pq->size >= pq->capacity){
		pq->capacity = 2*pq->capacity;
		pq->data = (node*) realloc(pq->data, 
                pq->capacity*sizeof(struct pqueue_node_t));
	}
	pq->data[pq->size-1] = create_node(key, ev);
  
	shiftup(pq, pq->size-1); 
} 
void pqueue_print(pqueue pq) {

    if (pqueue_empty(pq)) {
        printf("<empty queue>\n");
        return;
    }

    pqueue temp = pqueue_create();
    
    while (!pqueue_empty(pq)) {
        double key = pq->data[0]->key;
        event ev = pqueue_pop(pq);
        printf(" %.3lf: ", key);
        event_print(ev);
        pqueue_push(temp, key, ev);
    }
    while (!pqueue_empty(temp)) {
        double key = temp->data[0]->key;
        event ev = pqueue_pop(temp);
        pqueue_push(pq, key, ev);
    }
}
