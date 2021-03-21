#ifndef _UTILS_H_
#define _UTILS_H_

#include "tsp.h"

double max(double a, double b);
void print_error(const char *err);
double dist(int i, int j, instance inst);
int xpos(int i, int j, instance inst);
int xxpos(int i, int j, instance inst);
int upos(int i, instance inst);

#endif /* _UTILS_H_ */
