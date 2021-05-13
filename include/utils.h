#ifndef INCLUDE_UTILS_H_
#define INCLUDE_UTILS_H_

#include "../include/tsp.h"

/* cplex position helpers */
int xpos(int i, int j, int nnodes);
int xxpos(int i, int j, int nnodes);
int upos(int i, int nnodes);
int ypos(int i, int j, int nnodes);
edge xpos_inverse(int pos, int nnodes);

/* compute distances and zstar from tour */
double dist(int i, int j, instance inst);
double compute_distmatrix(instance inst);
double compute_zstar(instance inst, solution sol);

/* graphs utils */
int reachable(int* link, int i, int j);
int visitable(int* link, int nnodes);
int* edges_tosucc(edge* edges, int nnodes);

/* compare functions */
int wedgecmp(const void* a, const void* b);
int nodelexcmp(const void* a, const void* b);

/* computational geometry helpers */
double cross(node a, node b);
int ccw(node a, node b, node c);

/* quick string helpers */
char* model_type_tostring(enum model_types model_type);
char* model_folder_tostring(enum model_folders folder);

/* wall clock trackers */
int64_t stopwatch(struct timespec* s, struct timespec* e);
int64_t stopwatch_n(struct timespec* s, struct timespec* e);

/* classic utils */
double max(double a, double b);
double min(double a, double b);
int maxi(int a, int b);
int mini(int a, int b);
void swap(int* x, int* y);
void print_error(const char* err, ...);
char** list_files(enum model_folders folder, int* nmodels);

#endif  // INCLUDE_UTILS_H_
