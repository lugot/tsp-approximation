#ifndef INCLUDE_ADJLIST_H_
#define INCLUDE_ADJLIST_H_

typedef struct pair_t {
    int a, b;
} * pair;

typedef struct adjlist_t {
    pair* pairs;
    int npairs;
} * adjlist;

adjlist adjlist_create(int N);
void adjlist_free(adjlist l);
void adjlist_add_arc(adjlist l, int i, int j);
void adjlist_print(adjlist l);
pair* adjlist_loose_ends(adjlist l, int* npairs);


#endif  // INCLUDE_ADJLIST_H_
