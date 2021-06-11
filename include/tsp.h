#ifndef INCLUDE_TSP_H_
#define INCLUDE_TSP_H_

#include <cplex.h>

struct solution_t;

enum instance_types { TSP, TOUR, UNHANDLED_INSTANCE_TYPE };
enum costs { INTEGER, REAL };
typedef struct cplex_params_t {
    int randomseed;
    int num_threads;
    double timelimit;
    int available_memory;
    enum costs cost;
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
    char* model_name;
    char* model_comment;
    enum model_folders model_folder;
    enum instance_types instance_type;

    /* parameters */
    cplex_params params;

    /* data */
    enum weight_types weight_type;
    int nnodes;
    int ncols;
    node* nodes;
    double** adjmatrix;

    /* solutions */
    int nsols;
    struct solution_t** sols;
} * instance;

typedef struct doit_fn_input_t {
    int nedges;
    CPXCALLBACKCONTEXTptr context;
} * doit_fn_input;

enum model_types {
    OPTIMAL_TOUR,      // 0
    NOSEC,             // 1
    MTZ_STATIC,        // 2
    MTZ_LAZY,          // 3
    MTZ_INDICATOR,     // 4 -> MTZ: 28
    GGLIT_STATIC,      // 5
    GGLECT_STATIC,     // 6
    GGLIT_LAZY,        // 7
    BENDERS,           // 8
    BENDERS_CALLBACK,  // 9
    HARD_FIXING,       // 10
    SOFT_FIXING,       // 11
    MST,               // 12
    GREEDY,            // 13
    GRASP,             // 14
    EXTRA_MILEAGE,     // 15
    VNS,               // 16
    TABU_SEACH         // 17
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
void plot_profiler(instance* insts, int ninstances);

#endif  // INCLUDE_TSP_H_
