#ifndef INCLUDE_MODELS_BENDERS_H_
#define INCLUDE_MODELS_BENDERS_H_

#include <cplex.h>

#include "../../include/adjlist.h"
#include "../../include/tsp.h"

void perform_BENDERS(CPXENVptr env, CPXLPptr lp, instance inst, double* xstar,
                     struct timespec s, struct timespec e, int twophase);
int add_BENDERS_sec(CPXENVptr env, CPXLPptr lp, adjlist l);
int CPXPUBLIC add_BENDERS_sec_callback_driver(CPXCALLBACKCONTEXTptr context,
                                              CPXLONG contextid,
                                              void* userhandle);

#endif  // INCLUDE_MODELS_BENDERS_H_
