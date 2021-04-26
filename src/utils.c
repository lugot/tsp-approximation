#include "../include/utils.h"

#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include "../include/tsp.h"

double geodist(size_t i, size_t j, instance inst);
double l2dist(size_t i, size_t j, instance inst);

int xpos(int i, int j, int nnodes) {
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
    if (i > j) return xpos(j, i, nnodes);

    /*     area of triangle defined by row i */
    /*                       [^^^^^^^^^^^^^] */
    return i * nnodes + j - ((i + 1) * (i + 2)) / 2;
}
int xxpos(int i, int j, int nnodes) {
    /*
     * CPLEX variables representation:
     * in CPLEX variables are x_0, x_1, ... but I have i,j-indexed variables:
     * map needed
     *
     * for n=4 variables
     *   j 1 2 3 4
     * i
     * 1   0 1 2 3
     * 2   4 5 6 7
     * 3   8 9 10 11
     * 4   12 13 14 15
     */

    return i * nnodes + j;
}
int upos(int i, int nnodes) { return (nnodes * nnodes) + i - 1; }
int ypos(int i, int j, int nnodes) {
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

    int n = nnodes * nnodes + i * (nnodes - 1) + j;
    if (i < j) n--;

    return n;
}

double compute_zstar(instance inst, solution sol) {
    assert(inst->adjmatrix != NULL && "distances not computed yet");
    int nedges = sol->nedges;

    double zstar = 0.0;
    for (int i = 0; i < nedges; i++) {
        edge e = sol->edges[i];

        zstar += inst->adjmatrix[maxi(e.i, e.j)][mini(e.i, e.j)];
    }

    sol->zstar = zstar;
    return zstar;
}

/* distance functions helpers */
double geodist(size_t i, size_t j, instance inst) {
    // TODO(lugot): NOT CHECKED FOR FLOATING POINT SAFETY
    double RRR = 637.388;
    double q1 = cos(inst->nodes[i].x - inst->nodes[j].x);
    double q2 = cos(inst->nodes[i].y - inst->nodes[j].y);
    double q3 = cos(inst->nodes[i].y + inst->nodes[j].y);

    double distance =
        RRR * acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0;

    if (inst->params->cost == REAL) return distance;
    return 0.0 + (int)distance;
}

double l2dist(size_t i, size_t j, instance inst) {
    double dx = inst->nodes[i].x - inst->nodes[j].x;
    double dy = inst->nodes[i].y - inst->nodes[j].y;

    if (inst->params->cost == REAL) return sqrt(dx * dx + dy * dy);
    return 0.0 + (int)(sqrt(dx * dx + dy * dy) + 0.499999999);
}

double compute_dist(instance inst) {
    if (inst->adjmatrix != NULL) return 0.0;

    struct timeval start, end;
    gettimeofday(&start, NULL);

    inst->adjmatrix = (double**)calloc(inst->nnodes, sizeof(double*));

    for (int i = 0; i < inst->nnodes; i++) {
        /* allocate the lower row only */
        inst->adjmatrix[i] = (double*)calloc(i + 1, sizeof(double));

        for (int j = 0; j <= i; j++) {
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
    int64_t seconds = (end.tv_sec - start.tv_sec);
    int64_t micros = ((seconds * 1000000) + end.tv_usec) - (start.tv_usec);
    return micros / 1000.0;
}

double dist(int i, int j, instance inst) {
    assert(inst->adjmatrix != NULL && "distances not computed yet");

    return i > j ? inst->adjmatrix[i][j] : inst->adjmatrix[j][i];
}

/* check if node j is reachable in the list by node i */
int reachable(int* link, int i, int j) {
    int next = i;
    do {
        next = link[next];
        if (next == j) return 1;
    } while (link[next] != i);

    return 0;
}

/* check if all the nodes are visitated in a single tour*/
int visitable(int* link, int nnodes) {
    int visits = 0;

    int next = 0;
    do {
        next = link[next];
        visits++;
    } while (link[next] != 0);

    return visits == nnodes - 1;
}

/* print helpers */
void print_error(const char* err, ...) {
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
    int bufsize = 100;
    char* ans = (char*)calloc(bufsize, sizeof(char));

    switch (model_type) {
        case NOSEC:
            snprintf(ans, bufsize, "symmetric");
            break;
        case OPTIMAL_TOUR:
            snprintf(ans, bufsize, "optimal_tour");
            break;
        case MTZ_STATIC:
            snprintf(ans, bufsize, "mtz_static");
            break;
        case MTZ_LAZY:
            snprintf(ans, bufsize, "mtz_lazy");
            break;
        case GGLIT_STATIC:
            snprintf(ans, bufsize, "gglit_static");
            break;
        case GGLECT_STATIC:
            snprintf(ans, bufsize, "gglect_static");
            break;
        case GGLIT_LAZY:
            snprintf(ans, bufsize, "gglit_lazy");
            break;
        case BENDERS:
            snprintf(ans, bufsize, "benders");
            break;
        case BENDERS_CALLBACK:
            snprintf(ans, bufsize, "benders_callback");
            break;
        case HARD_FIXING:
            snprintf(ans, bufsize, "hard_fixing");
            break;
    }

    return ans;
}

/* random utils */
double max(double a, double b) { return a > b ? a : b; }
double min(double a, double b) { return a < b ? a : b; }
int maxi(int a, int b) { return a > b ? a : b; }
int mini(int a, int b) { return a < b ? a : b; }

void swap(int* x, int* y) {
    int temp;

    temp = *y;
    *y = *x;
    *x = temp;
}
