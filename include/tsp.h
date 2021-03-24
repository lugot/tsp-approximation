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
enum model_types {
	OPTIMAL_TOUR,
	SYMMETRIC,
	ASYMMETRIC_MTZ
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

    /* parameters */
    enum instance_types instance_type;
	int randomseed;
	int num_threads;
	double timelimit;
	int available_memory;
    enum costs costs_type;

	/* data */
	enum weight_types weight_type;
	int num_nodes;
	node* nodes;
	double** adjmatrix;

	/* solutions */
	int num_solutions;
	struct solution_t** sols;
} *instance;


typedef struct solution_t {
	struct instance_t* inst;
	enum model_types model_type;

	double zstar;
	int num_edges;
	edge* edges;
	int* parent;

	double distance_time;
	double build_time;
	double solve_time;
} *solution;



instance create_tsp_instance();
instance duplicate_instance(instance inst);
void free_instance();
void add_solution(instance inst, solution sol);

double build_tsp_model(instance inst, CPXENVptr env, CPXLPptr lp, enum model_types type);
void add_MTZ_subtour(instance inst, CPXENVptr env, CPXLPptr lp);

#endif   /* _TSP_H_ */
