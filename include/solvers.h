#ifndef _SOLVERS_H_
#define _SOLVERS_H_

#include "tsp.h"

solution TSPopt(instance inst, enum model_types model_type);
void save_results(instance* insts, int num_instances);

void retreive_symmetric_solution(double* xstar, solution sol);
void retreive_asymmetric_solution(double* xstar, solution sol);

#endif /* _SOLVERS_H_ */
