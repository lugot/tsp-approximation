#ifndef INCLUDE_UTILS_H_
#define INCLUDE_UTILS_H_

#include "../include/tsp.h"

int xpos(int i, int j, int nnodes);
int xxpos(int i, int j, int nnodes);
int upos(int i, int nnodes);
int ypos(int i, int j, int nnodes);
edge xpos_inverse(int pos, int nnodes);

double compute_zstar(instance inst, solution sol);

double dist(int i, int j, instance inst);
double l2dist(size_t i, size_t j, instance inst);
double compute_dist(instance inst);

int reachable(int* link, int i, int j);
int visitable(int* link, int nnodes);

int wedgecmp(const void* a, const void* b);
int nodelexcmp(const void* a, const void* b);

double cross(node a, node b);
int ccw(node a, node b, node c); 

char* model_type_tostring(enum model_types model_type);
char* model_folder_tostring(enum model_folders folder);

void print_error(const char* err, ...);

double max(double a, double b);
double min(double a, double b);
int maxi(int a, int b);
int mini(int a, int b);
void swap(int* x, int* y);

int64_t stopwatch(struct timespec* s, struct timespec* e);
int64_t stopwatch_n(struct timespec* s, struct timespec* e);

char** list_files(enum model_folders folder, int* nmodels);

#endif  // INCLUDE_UTILS_H_
