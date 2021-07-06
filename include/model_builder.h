#ifndef INCLUDE_MODEL_BUILDER_H_
#define INCLUDE_MODEL_BUILDER_H_

#include <cplex.h>

#include "../include/tsp.h"
#include "../include/adjlist.h"

double build_tsp_model(CPXENVptr env, CPXLPptr lp, instance inst,
                       enum model_types model_type);

#endif  // INCLUDE_MODEL_BUILDER_H_
