#ifndef INCLUDE_REFINEMENTS_H_
#define INCLUDE_REFINEMENTS_H_

#include "../include/tsp.h"

solution TSPtwoopt_multistart(instance inst);
solution TSPthreeopt_multistart(instance inst);

double twoopt_delta(instance inst, int* succ, int i, int j);

double twoopt_pick(instance inst, int* succ, int* a, int* b);
double twoopt_tabu_pick(instance inst, int* succ, int* tabu_nodes, int tenure,
                        int k, int* a, int* b);
void twoopt_move(int* succ, int nnodes, int a, int b);
double threeopt_pick(instance inst, int* succ, int* a, int* b, int* c,
                     struct timespec* s, struct timespec* e);
void threeopt_move(int* succ, int nnodes, int a, int b, int c, instance inst);

double twoopt_refinement_notimelim(instance inst, int* succ, int nnodes);
double twoopt_refinement(instance inst, int* succ, int nnodes,
                         struct timespec* s, struct timespec* e);
double threeopt_refinement(instance inst, int* succ, int nnodes,
                         struct timespec* s, struct timespec* e);
void kick(int* succ, int nnodes, int strength);

#endif  // INCLUDE_REFINEMENTS_H_
