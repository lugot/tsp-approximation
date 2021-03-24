#include "utils.h"
#include "tsp.h"
#include "math.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <sys/time.h>


double geodist(size_t i, size_t j, instance inst);
double l2dist(size_t i, size_t j, instance inst);

int xpos(int i, int j, instance inst) {
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
int upos(int i, instance inst) {
	return (inst->num_nodes * inst->num_nodes) + i-1;
}
int ypos(int i, int j, instance inst) {

	/*
	* for n=4 variables
	 *   j 1 2 3 4
	 * i
	 * 1   x 0 1 2
	 * 2   3 x 4 5
	 * 3   6 7 x 8
	 * 4   9 10 11 x
	 */

	if (i==j) print_error("variable y does not exist for same i and j");

	int n = inst->num_nodes * inst->num_nodes + i*(inst->num_nodes-1) + j;
	if (i < j) n--;

	return n;
}


double zstar(instance inst, solution sol) {
	if (inst->adjmatrix == NULL) {
		print_error("distances not computed yet");
	}

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

	if (inst->costs_type == REAL) return distance;
	return 0.0 + (int) distance;
}

double l2dist(size_t i, size_t j, instance inst) {
	double dx = inst->nodes[i].x - inst->nodes[j].x;
	double dy = inst->nodes[i].y - inst->nodes[j].y;

	if (inst->costs_type == REAL) return sqrt(dx*dx+dy*dy);
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
					print_error("unhandeled distance function");
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
	if (inst->adjmatrix == NULL) {
		print_error("distances not computed yet");
	}

	return i > j ? inst->adjmatrix[i][j] : inst->adjmatrix[j][i];
}



/* print helpers */

void print_error(const char *err, ...) {
    va_list args;
    va_start(args, err);

	printf("\n\n--- ERROR ---\n");
    vprintf(err, args);

    va_end( args );
	fflush(NULL);

	exit(EXIT_FAILURE);
}

void print_instance(instance inst, int print_data) {
	printf("\nmodel: %s\n", inst->model_name);
	printf("comment: %s\n", inst->model_comment);

	printf("parameters:\n");
	printf("- model type: ");
	switch (inst->instance_type) {
		case TSP:
			printf("TSP\n");
			break;
		case TOUR:
			printf("(solution only)\n");
			break;
		default:
			printf("\n");
			break;
	}
	printf("- random seed: %d\n", inst->randomseed);
	printf("- number of threads: %d\n", inst->num_threads);
	printf("- time limit: %lf\n", inst->timelimit);
	printf("- available memory: %d MB\n", inst->available_memory);
	printf("- costs type: ");
	switch (inst->costs_type) {
		case REAL:
			printf("real\n");
			break;
		case INTEGER:
			printf("integer\n");
			break;
	}

	printf("input data:\n");
	printf("- weight type: ");
	switch (inst->weight_type) {
		case ATT:
			printf("ATT\n");
			break;
		case EUC_2D:
			printf("EUC_2D\n");
			break;
		case GEO:
			printf("GEO\n");
			break;
		case EXPLICIT:
			printf("(matrix explicit)\n");
			break;
		default:
			break;
	}
	printf("- number of nodes: %d\n", inst->num_nodes);
	if (print_data) {
		printf("- weights:\n");
		if (inst->adjmatrix != NULL) {
			/* compute numof figures for spacing */
			int column_width = 0;
			for (int i=0; i<inst->num_nodes; i++) for (int j=0; j<=i; j++) {
				column_width = max(column_width, 1 + (int) log10(inst->adjmatrix[i][j]));
			}
			column_width = 2 + max(column_width, 1 + (int) log10(inst->num_nodes));

			/* figures.ab so max 1 decimal figures */
			char* buffer = (char*) calloc(column_width, sizeof(char));

			for (int i=0; i<inst->num_nodes; i++) {
				sprintf(buffer, "%d | ", i+1);
				printf("%*s", column_width+1, buffer);

				for (int j=0; j<=i; j++) {
					sprintf(buffer, "%.1lf ", inst->adjmatrix[i][j]);
					printf("%*s", column_width+1, buffer);
				}
				printf("\n");
			}
			free(buffer);
		}
		printf("- nodes:\n");
		if (inst->nodes != NULL) {
			int column_width = 0;
			for (int i=0; i<inst->num_nodes; i++) {
				column_width = max(column_width, 1 + (int) log10(inst->nodes[i].x));
				column_width = max(column_width, 1 + (int) log10(inst->nodes[i].y));
			}
			column_width = 2 + max(column_width, 1 + (int) log10(inst->num_nodes));

			/* figures.ab so max 1 decimal figures */
			char* buffer = (char*) calloc(column_width, sizeof(char));

			for (int i=0; i<inst->num_nodes; i++) {
				sprintf(buffer, "%d | ", i+1);
				printf("%*s", column_width+1, buffer);

				sprintf(buffer, "%.1lf ", inst->nodes[i].x);
				printf("%*s", column_width+1, buffer);

				sprintf(buffer, "%.1lf ", inst->nodes[i].y);
				printf("%*s\n", column_width+1, buffer);
			}
			free(buffer);
		}
	}

	printf("solutions:\n");
	printf("- num of solutions: %d\n", inst->num_solutions);
	for (int i=0; i<inst->num_solutions; i++) {
		printf("solution %d:\n", i+1);
		print_solution(inst->sols[i], print_data);
	}

	printf("--- ---\n\n");
}

void print_solution(solution sol, int print_data) {
	printf("- optimality: ");
	switch (sol->model_type) {
		case OPTIMAL_TOUR:
			printf("optimal tour from file\n");
			break;
		case SYMMETRIC:
			printf("symmetric, no subtour elimination\n");
			break;
		case ASYMMETRIC_MTZ:
			printf("asymmetric, mtz subtour elimination\n");
			break;
		case ASYMMETRIC_GG:
			printf("asymmetric, gg subtour elimination\n");
			break;
		default:
			printf("\n");
			break;
	}
	printf("- zstar: %lf\n", sol->zstar);
	printf("- num edges: %d\n", sol->num_edges);
	if (print_data) {
		printf("- edges:\n");
		if (sol->edges != NULL) {
			int column_width = 2 + (int) log10(sol->num_edges);
			char* buffer = (char*) calloc(column_width, sizeof(char));

			for (int i=0; i<sol->num_edges; i++) {
				sprintf(buffer, "%d | ", i+1);
				printf("%*s", column_width+2, buffer);

				sprintf(buffer, "%d ", sol->edges[i].i+1);
				printf("%*s", column_width, buffer);

				sprintf(buffer, "%d ", sol->edges[i].j+1);
				printf("%*s\n", column_width, buffer);
			}
			free(buffer);
		}
		printf("- parent:\n");
		if (sol->parent != NULL) {
			int column_width = 1 + (int) log10(sol->num_edges);
			char* buffer = (char*) calloc(column_width, sizeof(char));

			for (int i=0; i<sol->num_edges; i++) {
				sprintf(buffer, "%d", i+1);
				printf("%*s", column_width+1, buffer);
			}
			for (int i=0; i<sol->num_edges; i++) {
				sprintf(buffer, "%d", sol->parent[i]+1);
				printf("%*s", column_width+1, buffer);
			}
			printf("\n");
			free(buffer);
		}
	}
	printf("- distance time: %lf ms\n", sol->distance_time);
	printf("- build time: %lf ms\n", sol->build_time);
	printf("- solve time: %lf s\n", sol->solve_time / 1000.0);
}

void plot_solution_graphviz(solution sol) {
	plot_solutions_graphviz((solution[]) {sol}, 1);
}

void plot_solutions_graphviz(solution* sols, int num_sols) {
	double box_size = 20.0; // TODO: make proportional to the number of nodes
	double max_coord = 0.0;

	instance inst = sols[0]->inst;
	for (int i=0; i<inst->num_nodes && inst->instance_type == TSP; i++) {
		max_coord = max(max_coord, fabs(inst->nodes[i].x));
		max_coord = max(max_coord, fabs(inst->nodes[i].y));
	}

	char* filename;
	filename = (char*) calloc(100, sizeof(char));
	sprintf(filename, "../data/%s/%s.dot", inst->model_name, inst->model_name);

	FILE* fp;
	fp = fopen(filename, "w");

	fprintf(fp, "graph %s {\n", inst->model_name);
	fprintf(fp, "\tnode [shape=circle fillcolor=white]\n");
	for (int i=0; i<inst->num_nodes && inst->instance_type == TSP; i++) {
		double plot_posx = inst->nodes[i].x / max_coord * box_size;
		double plot_posy = inst->nodes[i].y / max_coord * box_size;

		fprintf(fp, "\t%d [ pos = \"%lf,%lf!\"]\n", i, plot_posx, plot_posy);
	}
	fprintf(fp, "\n");

	for (int u=0; u<num_sols; u++) {
		for (int k=0; k<sols[u]->num_edges; k++) {
			fprintf(fp, "\t%d -- %d", sols[u]->edges[k].i, sols[u]->edges[k].j);

			if (sols[u]->model_type == OPTIMAL_TOUR) fprintf(fp, " [color = red]");
			fprintf(fp, "\n");
		}
	}
	fprintf(fp, "}");

	fclose(fp);
	free(filename);
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
