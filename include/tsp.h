#ifndef _TSP_H_
#define _TSP_H_

enum tsp_model_type {NULL_MODEL, SYMMETRIC_TSP, TOUR};
enum tsp_costs_type {INTEGER, DOUBLE};

typedef struct tsp_instance_t {

    /* input data */
	int num_nodes;
	//double *demand;
	double* xcoord;
	double* ycoord;
	//int depot;
	//double capacity;
	//int nveh;
    char* model_name;
    char* model_comment;

    /* parameters */
    enum tsp_model_type model_type;
	char* input_file;
	//char* node_file;
	//int old_benders;
	int randomseed;
	int num_threads;
	double timelimit;
	int available_memory;
	//int max_nodes;
	//double cutoff;
    enum tsp_costs_type costs_type;

    /* global data */
	//double tstart;
	//double zbest;							// best sol. available
	//double tbest;							// time for the best sol. available
	//double* best_sol;						// best sol. available
	//double	best_lb;						// best lower bound available
	//double* load_min;						// minimum load when leaving a node
	//double* load_max;						// maximum load when leaving a node

    /* model */
	//int xstart;
	//int qstart;
	//int bigqstart;
	//int sstart;
	//int bigsstart;
	//int ystart;
	//int fstart;
	//int zstart;
} *tsp_instance;


/* public interface */
tsp_instance create_tsp_instance();
void print_instance(tsp_instance instance);

//inline
//inline int imax(int i1, int i2) { return ( i1 > i2 ) ? i1 : i2; }
//inline double dmin(double d1, double d2) { return ( d1 < d2 ) ? d1 : d2; }
//inline double dmax(double d1, double d2) { return ( d1 > d2 ) ? d1 : d2; }

#endif   /* _TSP_H_ */
