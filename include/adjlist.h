#ifndef INCLUDE_ADJLIST_H_
#define INCLUDE_ADJLIST_H_

typedef struct adjlist_t {
    int** neigh;
    int N;
    int* visited;
    int i;
    int bit;
} * adjlist;

adjlist adjlist_create(int N);
void adjlist_add_arc(adjlist l, int i, int j);
int adjlist_get_arc(adjlist l, int* i, int* j);
int adjlist_loose_ends(adjlist l, int* end1, int* end2);
void adjlist_print(adjlist l);
void adjlist_reset(adjlist l);
void adjlist_free(adjlist l);

#endif  // INCLUDE_ADJLIST_H_
