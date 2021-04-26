#ifndef INCLUDE_UTILS_H_
#define INCLUDE_UTILS_H_

#include "../include/tsp.h"

int xpos(int i, int j, int nnodes);
int xxpos(int i, int j, int nnodes);
int upos(int i, int nnodes);
int ypos(int i, int j, int nnodes);

double compute_zstar(instance inst, solution sol);

double dist(int i, int j, instance inst);
double compute_dist(instance inst);

int reachable(int* link, int i, int j);
int visitable(int* link, int nnodes);

char* model_type_tostring(enum model_types model_type);

void print_error(const char* err, ...);

double max(double a, double b);
double min(double a, double b);
int maxi(int a, int b);
int mini(int a, int b);
void swap(int* x, int* y);

#endif  // INCLUDE_UTILS_H_
