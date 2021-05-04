#ifndef INCLUDE_ADJLIST_H_
#define INCLUDE_ADJLIST_H_

typedef struct adjlist_t {
    int (*neigh)[2];
    int N;
    int* visited;
    int i;
    int bit;
} * adjlist;

adjlist adjlist_create(int N);
void adjlist_add_edge(adjlist l, int i, int j);
int adjlist_get_edge(adjlist l, int* i, int* j);
int adjlist_get_loose_ends(adjlist l, int* end1, int* end2);
int* adjlist_get_subtour(adjlist l, int* subsize);
int adjlist_single_tour(adjlist l);
void adjlist_print(adjlist l);
void adjlist_reset(adjlist l);
void adjlist_hard_reset(adjlist l);
void adjlist_free(adjlist l);

#endif  // INCLUDE_ADJLIST_H_
