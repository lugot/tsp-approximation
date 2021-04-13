#ifndef _MODEL_BUILDER_H_
#define _MODEL_BUILDER_H_

#include "tsp.h"
#include <cplex.h>

double build_tsp_model(CPXENVptr env, CPXLPptr lp, instance inst, enum model_types model_type);
void add_BENDERS_sec(CPXENVptr env, CPXLPptr lp, solution sol);

#endif /* _MODEL_BUILDER_H_ */
