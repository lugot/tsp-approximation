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
	SYMMETRIC_SUBTOUR,
	ASYMMETRIC_MTZ
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
	int i, j;
} edge;


typedef struct instance_t {
	/* model infos */
    char* model_name;
    char* model_comment;

    /* parameters */
    enum types type;
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

	enum optimalities optimality;
	enum tsp_types tsp_type;

	double zstar;

	int num_edges;
	edge* edges;
	int* parent;

} *solution;



instance create_tsp_instance();
instance duplicate_instance(instance inst);
void free_instance();
void add_solution(instance inst, solution sol);

void build_tsp_model(instance inst, CPXENVptr env, CPXLPptr lp, enum optimalities opt);
void add_MTZ_subtour(instance inst, CPXENVptr env, CPXLPptr lp);

#endif   /* _TSP_H_ */
