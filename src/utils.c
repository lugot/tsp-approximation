#include "utils.h"
#include "tsp.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <assert.h>
#include <sys/time.h>


double geodist(size_t i, size_t j, instance inst);
double l2dist(size_t i, size_t j, instance inst);

int xpos(int i, int j, int num_nodes) {
	/*
 	 * CPLEX variables representation:
 	 * in CPLEX variables are x_0, x_1, ... but I have i,j-indexed variables:
 	 * map needed
 	 *
 	 * for n=4 variables
 	 *   j 1 2 3 4
 	 * i
 	 * 1   X 0 1 2
 	 * 2   X X 3 4
 	 * 3   X X X X
 	 * 4   X X X X
 	 */

	if (i == j) print_error("i == j in xpos");
	if (i > j) return xpos(j, i, num_nodes);

	/*     area of triangle defined by row i */
	/*                       [^^^^^^^^^^^^^] */
	return i*num_nodes + j - ((i+1)*(i+2))/2;
}
int xxpos(int i, int j, int num_nodes) {
	/*
	 * CPLEX variables representation:
	 * in CPLEX variables are x_0, x_1, ... but I have i,j-indexed variables: map needed
	 *
	 * for n=4 variables
	 *   j 1 2 3 4
	 * i
	 * 1   0 1 2 3
	 * 2   4 5 6 7
	 * 3   8 9 10 11
	 * 4   12 13 14 15
	 */

	return i*num_nodes + j;
}
int upos(int i, int num_nodes) {
	return (num_nodes * num_nodes) + i-1;
}
int ypos(int i, int j, int num_nodes) {
	/*
	* for n=4 variables
	 *   j 1 2 3 4
	 * i
	 * 1   x 0 1 2
	 * 2   3 x 4 5
	 * 3   6 7 x 8
	 * 4   9 10 11 x
	 */

	if (i == j) print_error("variable y does not exist for same i and j");

	int n = num_nodes * num_nodes + i*(num_nodes-1) + j;
	if (i < j) n--;

	return n;
}


double zstar(instance inst, solution sol) {
	assert(inst->adjmatrix != NULL && "distances not computed yet");

	double zstar = 0.0;
	for (int i=0; i<sol->num_edges; i++) {
		edge e = sol->edges[i];

		zstar += inst->adjmatrix[maxi(e.i, e.j)][mini(e.i, e.j)];
	}

	sol->zstar = zstar;
	return zstar;
}


/* distance functions helpers */
double geodist(size_t i, size_t j, instance inst) {
	// TODO: NOT CHECKED FOR FLOATING POINT SAFETY
	double RRR = 637.388;
	double q1 = cos(inst->nodes[i].x - inst->nodes[j].x);
	double q2 = cos(inst->nodes[i].y - inst->nodes[j].y);
	double q3 = cos(inst->nodes[i].y + inst->nodes[j].y);

	double distance = RRR * acos(0.5*((1.0+q1)*q2 - (1.0-q1)*q3)) + 1.0;

	if (inst->params->cost == REAL) return distance;
	return 0.0 + (int) distance;
}

double l2dist(size_t i, size_t j, instance inst) {
	double dx = inst->nodes[i].x - inst->nodes[j].x;
	double dy = inst->nodes[i].y - inst->nodes[j].y;

	if (inst->params->cost == REAL) return sqrt(dx*dx+dy*dy);
	return 0.0 + (int) (sqrt(dx*dx+dy*dy) + 0.499999999);
}

double compute_dist(instance inst) {
	if (inst->adjmatrix != NULL) return 0.0;

	struct timeval start, end;
    gettimeofday(&start, NULL);

	inst->adjmatrix = (double**) calloc(inst->num_nodes, sizeof(double*));

	for (int i=0; i<inst->num_nodes; i++) {
		/* allocate the lower row only */
		inst->adjmatrix[i] = (double*) calloc(i+1, sizeof(double));

		for (int j=0; j<=i; j++) {
			switch (inst->weight_type) {
				case ATT:
				case EUC_2D:
					inst->adjmatrix[i][j] = l2dist(i, j, inst);
					break;

				case GEO:
					inst->adjmatrix[i][j] = geodist(i, j, inst);
					break;

				default:
					inst->adjmatrix[i][j] = 0.0;
					assert(1 && "unhandeled distance function");
					break;
			}
		}
	}

    gettimeofday(&end, NULL);
    long seconds = (end.tv_sec - start.tv_sec);
    long micros = ((seconds * 1000000) + end.tv_usec) - (start.tv_usec);
	return micros / 1000.0;
}

double dist(int i, int j, instance inst) {
	assert(inst->adjmatrix != NULL && "distances not computed yet");

	return i > j ? inst->adjmatrix[i][j] : inst->adjmatrix[j][i];
}

/* check if node j is reachable in the list by node i */
int reachable(solution sol, int i, int j) {
	int next = i;
	do {
		next = sol->parent[next];
		if (next == j) return 1;
	} while (sol->parent[next] != i);

	return 0;
}

/* check if all the nodes are visitated in a single tour*/
int visitable(solution sol) {
	int visits = 0;

	int next = 0;
	do {
		next = sol->parent[next];
		visits++;
	} while (sol->parent[next] != 0);

	return visits == sol->num_edges-1;
}


/* print helpers */
void print_error(const char *err, ...) {
    va_list args;
    va_start(args, err);

	printf("\n\n--- ERROR ---\n");
    vprintf(err, args);
	printf("\n");

    va_end(args);
	fflush(NULL);

	exit(EXIT_FAILURE);
}

char* model_type_tostring(enum model_types model_type) {
	char* ans = (char*) calloc(20, sizeof(char));

	switch (model_type) {
		case SYMMETRIC:
			strcpy(ans, "symmetric");
			break;

		case OPTIMAL_TOUR:
			strcpy(ans, "optimal_tour");
			break;

		case ASYMMETRIC_MTZ:
			strcpy(ans, "mtz");
			break;

		case ASYMMETRIC_GG:
			strcpy(ans, "gg");
			break;

		case SYMMETRIC_BENDERS:
			strcpy(ans, "benders");
			break;

		case SYMMETRIC_BENDERS_CALLBACK: ;
			strcpy(ans, "benders_callback");
			break;
	}

	return ans;
}


/* random utils */
double max(double a, double b) {
	return a > b ? a : b;
}
double min(double a, double b) {
	return a < b ? a : b;
}
int maxi(int a, int b) {
	return a > b ? a : b;
}
int mini(int a, int b) {
	return a < b ? a : b;
}

void swap(int* x, int* y) {
   int temp;

   temp = *y;
   *y = *x;
   *x = temp;
}
