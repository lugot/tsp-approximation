#include "utils.h"
#include "tsp.h"
#include "math.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>


double geodist(size_t i, size_t j, instance inst);
double l2dist(size_t i, size_t j, instance inst);
void compute_dist(instance inst);

void print_error(const char *err, ...) {
    va_list args;
    va_start(args, err);

	printf("\n\n--- ERROR ---\n");
    vprintf(err, args);

    va_end( args );
	fflush(NULL);

	exit(EXIT_FAILURE);
}

double geodist(size_t i, size_t j, instance inst) {

	double dx = inst->nodes[i].x - inst->nodes[j].x;
	double dy = inst->nodes[i].y - inst->nodes[j].y;

	if (inst->costs_type == REAL) return sqrt(dx*dx+dy*dy);
	return 0.0 + (int) (sqrt(dx*dx+dy*dy) + 0.499999999);
}
double l2dist(size_t i, size_t j, instance inst) {

	double dx = inst->nodes[i].x - inst->nodes[j].x;
	double dy = inst->nodes[i].y - inst->nodes[j].y;

	if (inst->costs_type == REAL) return sqrt(dx*dx+dy*dy);
	return 0.0 + (int) (sqrt(dx*dx+dy*dy) + 0.499999999);
}
void compute_dist(instance inst) {

	inst->adjmatrix = (double**) calloc(inst->num_nodes, sizeof(double*));

	for (int i=0; i<inst->num_nodes; i++) {
		/* allocate the lower row only */
		inst->adjmatrix[i] = (double*) calloc(i+1, sizeof(double));

		for (int j=0; j<=i; j++) {
			double distance;

			switch (inst->weight_type) {
				case ATT:
				case EUC_2D:
					distance = l2dist(i, j, inst);
					break;

				case GEO:
					distance = geodist(i, j, inst);
					break;

				default:
					distance = 0.0;
					print_error("unhandeled distance function");
					break;
			}

			inst->adjmatrix[i][j] = distance;
		}
	}

}
double dist(int i, int j, instance inst) {
	if (inst->adjmatrix == NULL) compute_dist(inst);

	return i > j ? inst->adjmatrix[i][j] : inst->adjmatrix[j][i];
}

/*int xpos(int i, int j, instance inst) {*/
	/*
	 * CPLEX variables representation:
	 * in CPLEX variables are x_0, x_1, ... but I have i,j-indexed variables: map needed
	 *
	 * for n=4 variables
	 *   j 1 2 3 4
	 * i
	 * 1   X 0 1 2
	 * 2   X X 3 4
	 * 3   X X X X
	 * 4   X X X X
	 */

	/*if (i == j) print_error("i == j in xpos");*/
	/*if (i > j) return xpos(j, i, inst);*/

	/*[>           area of triangle defined by row i <]*/
	/*[>                             [^^^^^^^^^^^^^] <]*/
	/*return i*inst->num_nodes + j - ((i+1)*(i+2))/2;*/
/*}*/


double max(double a, double b) {
	return a > b ? a : b;
}
double min(double a, double b) {
	return a < b ? a : b;
}
