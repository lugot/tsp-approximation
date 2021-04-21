#ifndef _MODEL_BUILDER_H_
#define _MODEL_BUILDER_H_

#include <cplex.h>

#include "tsp.h"

double build_tsp_model(CPXENVptr env, CPXLPptr lp, instance inst,
                       enum model_types model_type);
void add_BENDERS_sec(CPXENVptr env, CPXLPptr lp, solution sol);
int CPXPUBLIC add_BENDERS_sec_callback_driver(CPXCALLBACKCONTEXTptr context,
                                              CPXLONG contextid,
                                              void *userhandle);

#endif /* _MODEL_BUILDER_H_ */
