#ifndef INCLUDE_ADJLIST_H_
#define INCLUDE_ADJLIST_H_

typedef struct pair_t {
    int a, b;
} * pair;

typedef struct adjlist_t {
    pair* pairs;
    int npairs;
    int* visited;
    int i;
} * adjlist;

adjlist adjlist_create(int N);
void adjlist_free(adjlist l);
void adjlist_add_arc(adjlist l, int i, int j);
void adjlist_print(adjlist l);
int adjlist_loose_ends(adjlist l, int* end1, int* end2);


#endif  // INCLUDE_ADJLIST_H_
