#ifndef INCLUDE_REFINEMENTS_H_
#define INCLUDE_REFINEMENTS_H_

#include "../include/tsp.h"

double twoopt_refinement(instance inst, int* succ, int nnodes);
double threeopt_refinement(instance inst, int* succ, int nnodes);
void kick(int* succ, int nnodes, int strength);

#endif  // INCLUDE_REFINEMENTS_H_
