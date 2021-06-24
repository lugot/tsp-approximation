#ifndef INCLUDE_MODELS_FIXING_H_
#define INCLUDE_MODELS_FIXING_H_

#include <cplex.h>

#include "../include/tsp.h"

void perform_HARD_FIXING(CPXENVptr env, CPXLPptr lp, instance inst,
                         double* xstar, int perc);
void perform_SOFT_FIXING(CPXENVptr env, CPXLPptr lp, instance inst,
                         double* xstar);

#endif  // INCLUDE_MODELS_FIXING_H_
