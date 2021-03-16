#ifndef _TSP_H_
#define _TSP_H_

#include <cplex.h>

enum optimalities {
	OPTIMAL_TOUR,
	OPTIMAL,
	APPROXIMATED
};

typedef struct edge_t {
	size_t i, j;
} edge;

typedef struct solution_t {
	enum optimalities optimality;
	double zstar;
	size_t num_edges;
	edge* edges;
} *solution;


enum model_types {
	NULL_MODEL,
	SYMMETRIC_TSP
};
enum costs_types {
	INTEGER,
	REAL
};

typedef struct node_t {
	double x, y;
} node;

typedef struct instance_t {
    /* input data */
	size_t num_nodes;
	node* nodes;

	/* model infos */
    char* model_name;
    char* model_comment;

    /* parameters */
    enum model_types model_type;
	int randomseed;
	int num_threads;
	double timelimit;
	int available_memory;
    enum costs_types costs_type;
} *instance;


instance create_tsp_instance();
void build_tsp_model(instance inst, CPXENVptr env, CPXLPptr lp);
void print_instance(instance inst);
void print_solution(solution sol);
void plot_solutions_graphviz(instance inst, solution* sols, size_t num_sols);
void plot_solution_graphviz(instance inst, solution sol);

#endif   /* _TSP_H_ */
