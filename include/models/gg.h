#ifndef INCLUDE_MODELS_GG_H_
#define INCLUDE_MODELS_GG_H_

#include <cplex.h>

#include "../include/model_builder.h"

void add_GG_variables(CPXENVptr env, CPXLPptr lp, instance inst);
void add_GGlit_static_sec(CPXENVptr env, CPXLPptr lp, instance inst);
void add_GGlect_static_sec(CPXENVptr env, CPXLPptr lp, instance inst);
void add_GGlit_lazy_sec(CPXENVptr env, CPXLPptr lp, instance inst);
void add_GGlect_lazy_sec(CPXENVptr env, CPXLPptr lp, instance inst);

#endif  // INCLUDE_MODELS_GG_H_
