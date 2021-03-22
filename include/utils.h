#ifndef _UTILS_H_
#define _UTILS_H_

#include "tsp.h"

int xpos(int i, int j, instance inst);
int xxpos(int i, int j, instance inst);
int upos(int i, instance inst);

double dist(int i, int j, instance inst);

void print_error(const char *err, ...);
void print_instance(instance inst);
void print_solution(solution sol);
void plot_solution_graphviz(solution sol);
void plot_solutions_graphviz(solution* sols, int num_sols);

double max(double a, double b);
double min(double a, double b);

#endif /* _UTILS_H_ */
