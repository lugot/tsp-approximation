#ifndef INCLUDE_REFINEMENTS_H_
#define INCLUDE_REFINEMENTS_H_

#include "../include/tsp.h"

double twoopt_delta(instance inst, int* succ, int i, int j);

double twoopt_pick(instance inst, int* succ, int* a, int* b);
double twoopt_tabu_pick(instance inst, int* succ, int* tabu_nodes, int tenure,
                        int k, int* a, int* b);
void twoopt_move(int* succ, int nnodes, int a, int b);
double threeopt_pick(instance inst, int* succ, int* a, int* b, int* c);
void threeopt_move(int* succ, int nnodes, int a, int b, int c, instance inst);

double twoopt_refinement(instance inst, int* succ, int nnodes);
double threeopt_refinement(instance inst, int* succ, int nnodes);
void kick(int* succ, int nnodes, int strength);

#endif  // INCLUDE_REFINEMENTS_H_
