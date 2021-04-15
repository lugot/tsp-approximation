#ifndef _TSP_H_
#define _TSP_H_

#include <cplex.h>

struct solution_t;

enum instance_types {
	TSP,
	TOUR,
	UNHANDLED_INSTANCE_TYPE
};
enum costs {
	INTEGER,
	REAL
};
typedef struct cplex_params_t {
	int randomseed;
	int num_threads;
	double timelimit;
	int available_memory;
    enum costs cost;

} *cplex_params;


enum model_folders {
	TSPLIB,
	GENERATED
};
enum weight_types {
	ATT,
	EUC_2D,
	GEO,
	EXPLICIT,
	UNHANDLED_WEIGHT_TYPE
};
typedef struct node_t {
	double x, y;
} node;
typedef struct edge_t {
	int i, j;
} edge;
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
	int num_nodes;
	node* nodes;
	double** adjmatrix;

	/* solutions */
	int num_solutions;
	struct solution_t** sols;
}
*instance;


enum model_types {
	OPTIMAL_TOUR,
	SYMMETRIC,
	ASYMMETRIC_MTZ,
	ASYMMETRIC_GG,
	ASYMMETRIC_PROF_GG,
	SYMMETRIC_BENDERS,
	SYMMETRIC_BENDERS_CALLBACK
};
typedef struct solution_t {
	struct instance_t* inst;
	enum model_types model_type;

	double zstar;
	int num_edges;
	edge* edges;
	int* parent;

	double start;
	double end;
	double distance_time;
	double build_time;
	double solve_time;
} *solution;

/* instance manipulators */
instance create_empty_instance();
instance create_instance(cplex_params params);
void add_params(instance inst, cplex_params params);
instance generate_random_instance(int id, int num_nodes, int box_size);
instance* generate_random_instances(int num_instances, int num_nodes, int box_size);
instance clone_instance(instance inst);
void free_instance();
void add_solution(instance inst, solution sol);

/* printers */
void print_instance(instance inst, int print_data);
void print_cplex_params(cplex_params params);
void print_solution(solution sol, int print_data);
void plot_solution_graphviz(solution sol);
void plot_solutions_graphviz(solution* sols, int num_sols);

#endif /* _TSP_H_ */
