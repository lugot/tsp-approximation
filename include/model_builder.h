#ifndef INCLUDE_MODEL_BUILDER_H_
#define INCLUDE_MODEL_BUILDER_H_

#include <cplex.h>

#include "../include/tsp.h"
#include "../include/adjlist.h"

double build_tsp_model(CPXENVptr env, CPXLPptr lp, instance inst,
                       enum model_types model_type);
void perform_BENDERS(CPXENVptr env, CPXLPptr lp, instance inst, double* xstar,
                     struct timespec s, struct timespec e);
void add_BENDERS_sec(CPXENVptr env, CPXLPptr lp, adjlist l);
void perform_HARD_FIXING(CPXENVptr env, CPXLPptr lp, instance inst,
                         double* xstar, int perc);
void perform_SOFT_FIXING(CPXENVptr env, CPXLPptr lp, instance inst,
                         double* xstar);
int CPXPUBLIC add_BENDERS_sec_callback_driver(CPXCALLBACKCONTEXTptr context,
                                              CPXLONG contextid,
                                              void *userhandle);

#endif  // INCLUDE_MODEL_BUILDER_H_
