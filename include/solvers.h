#ifndef INCLUDE_SOLVERS_H_
#define INCLUDE_SOLVERS_H_

#include "../include/tsp.h"

solution TSPopt(instance inst, enum model_types model_type);

void get_symmsol(double* xstar, solution sol);
void get_asymmsol(double* xstar, solution sol);

#endif  // INCLUDE_SOLVERS_H_
