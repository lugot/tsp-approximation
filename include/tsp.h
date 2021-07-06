#ifndef INCLUDE_TSP_H_
#define INCLUDE_TSP_H_

#include <cplex.h>

#include "../include/tracker.h"

struct solution_t;

typedef struct cplex_params_t {
    int randomseed;
    int num_threads;
    double timelimit;
    int available_memory;
} * cplex_params;

enum model_folders { TSPLIB, GENERATED };
enum weight_types { ATT, EUC_2D, GEO, EXPLICIT, UNHANDLED_WEIGHT_TYPE };
typedef struct node_t {
    double x, y;
} node;
typedef struct edge_t {
    int i, j;
} edge;
typedef struct wedge_t {
    double w;
    int i, j;
} wedge;
typedef struct instance_t {
    /* model infos */
    char* instance_name;
    char* instance_comment;
    char* instance_folder;

    /* parameters */
    cplex_params params;

    /* data */
    enum weight_types weight_type;
    int nnodes;
    int ncols;
    node* nodes;

    /* solutions */
    double zbest;
    int nsols;
    struct solution_t** sols;
} * instance;


enum model_types {
    OPTIMAL_TOUR,
    NOSEC,
    MTZ_STATIC,
    MTZ_LAZY,
    MTZ_LAZY_DEG2,
    MTZ_LAZY_DEG3,
    MTZ_INDICATOR,
    GGLIT_STATIC,
    GGLECT_STATIC,
    GGLIT_LAZY,
    GGLECT_LAZY,
    GGLIT_STATIC_DEG2,
    BENDERS,
    BENDERS_TWOPHASES,
    BENDERS_CALLBACK,
    HARD_FIXING,
    SOFT_FIXING,
    MST,
    GREEDY,
    GRASP,
    EXTRA_MILEAGE,
    TWOOPT_MULTISTART,
    VNS_RANDOM,
    VNS_GREEDY,
    TABU_SEACH
};
typedef struct solution_t {
    struct instance_t* inst;
    enum model_types model_type;

    double zstar;
    int nedges;
    edge* edges;

    double start;
    double end;
    double distance_time;
    double build_time;
    double solve_time;

    tracker t;
} * solution;

/* cplex param manipulators */
cplex_params create_params();
void add_params(instance inst, cplex_params params);
void free_params(cplex_params params);

/* instance manipulators */
instance create_empty_instance();
instance create_instance(cplex_params params);
instance generate_random_instance(int id, int num_nodes);
instance* generate_random_instances(int num_instances, int num_nodes);
void save_instance(instance inst);
void free_instance();

/* solution manipulators */
solution create_solution(instance inst, enum model_types model_type,
                         int nedges);
void add_solution(instance inst, solution sol);
void free_solution(solution sol);

/* printers */
void print_instance(instance inst, int print_data);
void print_cplex_params(cplex_params params);
void print_solution(solution sol, int print_data);

/* plotters */
void plot_graphviz(solution sol, int* edgecolors, int version);
void plot_profiler(instance* insts, int ninstances, int plot_obj);
void plot_tracking(instance* insts, int ninstances, int dist);

#endif  // INCLUDE_TSP_H_
