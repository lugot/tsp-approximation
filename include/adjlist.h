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
void add_arc(adjlist l, int a, int b);


#endif  // INCLUDE_ADJLIST_H_
