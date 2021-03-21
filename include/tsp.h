#ifndef _TSP_H_
#define _TSP_H_

#include <cplex.h>

struct solution_t;

enum types {
	TSP,
	TOUR,
	UNMANAGED_TYPE
};
enum costs {
	INTEGER,
	REAL
};
enum optimalities {
	OPTIMAL_TOUR,
	OPTIMAL,
	APPROXIMATED
};
enum tsp_types {
	SYMMETRIC,
	ASYMMETRIC
};
enum weight_types {
	ATT,
	EUC_2D,
	GEO,
	EXPLICIT,
	UNMANAGED_WEIGHT_TYPE
};


typedef struct node_t {
	double x, y;
} node;
typedef struct edge_t {
	size_t i, j;
} edge;


typedef struct instance_t {
	/* model infos */
    char* model_name;
    char* model_comment;

    /* parameters */
    enum types model_type;
	int randomseed;
	int num_threads;
	double timelimit;
	int available_memory;
    enum costs costs_type;

	/* data */
	enum weight_types weight_type;
	size_t num_nodes;
	node* nodes;
	double** adjmatrix;

	/* solutions */
	size_t num_solutions;
	struct solution_t** sols;
} *instance;


typedef struct solution_t {
	struct instance_t* inst;

	enum optimalities optimality;
	enum tsp_types tsp_type;


	double zstar;

	size_t num_edges;
	edge* edges;
	size_t* parent;

} *solution;



instance create_tsp_instance();
void build_tsp_model(instance inst, CPXENVptr env, CPXLPptr lp);
void print_instance(instance inst);
void print_solution(solution sol);
void plot_solutions_graphviz(instance inst, solution* sols, size_t num_sols);
void plot_solution_graphviz(instance inst, solution sol);

#endif   /* _TSP_H_ */
