#include "utils.h"
#include "tsp.h"
#include "math.h"
#include <stdlib.h>
#include <stdio.h>

double max(double a, double b) {
	return a > b ? a : b;
}

void print_error(const char *err) {
	printf("\n\n ERROR: %s \n\n", err);
	fflush(NULL);

	exit(EXIT_FAILURE);
}

double dist(int i, int j, instance inst) {
	/*
	 *l2 distance between nodes coordinates
	 */

	double dx = inst->nodes[i].x - inst->nodes[j].x;
	double dy = inst->nodes[i].y - inst->nodes[j].y;

	if (inst->costs_type == REAL) return sqrt(dx*dx+dy*dy);
	return 0.0 + (int) (sqrt(dx*dx+dy*dy) + 0.499999999);
}

int xpos(int i, int j, instance inst) {
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

	if (i == j) print_error("i == j in xpos");
	if (i > j) return xpos(j, i, inst);

	/*           area of triangle defined by row i */
	/*                             [^^^^^^^^^^^^^] */
	return i*inst->num_nodes + j - ((i+1)*(i+2))/2;
}

int xxpos(int i, int j, instance inst) {
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

	return i*inst->num_nodes + j;
}
