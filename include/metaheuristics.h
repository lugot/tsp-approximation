#ifndef INCLUDE_METAHEURISTICS_H_
#define INCLUDE_METAHEURISTICS_H_

#include "../include/tsp.h"

solution TSPvns(instance inst, int* succ);
solution TSPtabusearch(instance inst, int* succ);

#endif  // INCLUDE_METAHEURISTICS_H_

