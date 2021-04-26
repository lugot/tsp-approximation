#ifndef INCLUDE_SOLVERS_H_
#define INCLUDE_SOLVERS_H_

#include "../include/tsp.h"

solution TSPopt(instance inst, enum model_types model_type);
void save_results(instance* insts, int num_instances);

void get_symmsol(double* xstar, int nedges, edge* edges, int* link);
void get_asymmsol(double* xstar, int nedges, edge* edges, int* link);

#endif  // INCLUDE_SOLVERS_H_
