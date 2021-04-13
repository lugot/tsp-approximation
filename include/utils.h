#ifndef _UTILS_H_
#define _UTILS_H_

#include "tsp.h"

int xpos(int i, int j, int num_nodes);
int xxpos(int i, int j, int num_nodes);
int upos(int i, int num_nodes);
int ypos(int i, int j, int num_nodes);

double zstar(instance inst, solution sol);

double dist(int i, int j, instance inst);
double compute_dist(instance inst);

int reachable(solution sol, int i, int j);
int visitable(solution sol);

void print_error(const char *err, ...);

double max(double a, double b);
double min(double a, double b);
int maxi(int a, int b);
int mini(int a, int b);
void swap(int* x, int* y);

#endif /* _UTILS_H_ */
