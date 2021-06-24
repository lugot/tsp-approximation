#ifndef INCLUDE_MODELS_MTZ_H_
#define INCLUDE_MODELS_MTZ_H_

#include <cplex.h>

#include "../include/tsp.h"

void add_MTZ_variables(CPXENVptr env, CPXLPptr lp, instance inst);
void add_MTZ_static_sec(CPXENVptr env, CPXLPptr lp, instance inst);
void add_MTZ_lazy_sec(CPXENVptr env, CPXLPptr lp, instance inst);
void add_MTZ_indicator_sec(CPXENVptr env, CPXLPptr lp, instance inst);

#endif  // INCLUDE_MODELS_MTZ_H_
