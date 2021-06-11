#include "../include/utils.h"

#include <assert.h>
#include <dirent.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

#include "../include/globals.h"
#include "../include/tsp.h"

double l2dist(size_t i, size_t j, instance inst);
double geodist(size_t i, size_t j, instance inst);

/* cplex position helpers */
int xpos(int i, int j, int nnodes) {
    /*
     * CPLEX variables representation:
     * in CPLEX variables are x_0, x_1, ... but I have i,j-indexed variables:
     * map needed
     *
     * for n=4 variables
     *   j 0 1 2 3
     * i
     * 0   X 0 1 2
     * 1   X X 3 4
     * 2   X X X 5
     * 3   X X X X
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
     *   j 0 1 2 3
     * i
     * 0   0 1 2 3
     * 1   4 5 6 7
     * 2   8 9 10 11
     * 3   12 13 14 15
     */

    return i * nnodes + j;
}
int upos(int i, int nnodes) { return (nnodes * nnodes) + i - 1; }
int ypos(int i, int j, int nnodes) {
    /*
     * for n=4 variables
     *   j 0 1 2 3
     * i
     * 0   x 0 1 2
     * 1   3 x 4 5
     * 2   6 7 x 8
     * 3   9 10 11 x
     */

    if (i == j) print_error("variable y does not exist for same i and j");

    int n = nnodes * nnodes + i * (nnodes - 1) + j;
    if (i < j) n--;

    return n;
}
edge xpos_inverse(int pos, int nnodes) {
    // TODO(lugot): PERFORMANCE
    int tosub = nnodes - 1;
    int temp = pos;

    edge e;
    e.i = 0;
    while (temp >= tosub) {
        temp -= tosub;
        tosub--;
        e.i++;
    }
    e.j = e.i + 1 + temp;

    return e;
}

/* compute distances and zstar from tour */
double l2dist(size_t i, size_t j, instance inst) {
    double dx = inst->nodes[i].x - inst->nodes[j].x;
    double dy = inst->nodes[i].y - inst->nodes[j].y;

    if (inst->params->cost == REAL) return sqrt(dx * dx + dy * dy);
    return 0.0 + (int)(sqrt(dx * dx + dy * dy) + 0.499999999);
}
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
double dist(int i, int j, instance inst) {
    if (inst->adjmatrix != NULL) {
        return i > j ? inst->adjmatrix[i][j] : inst->adjmatrix[j][i];
    }

    switch (inst->weight_type) {
        case ATT:
        case EUC_2D:
            return l2dist(i, j, inst);
            break;
        case GEO:
            return geodist(i, j, inst);
            break;
        case EXPLICIT:
        case UNHANDLED_WEIGHT_TYPE:
            assert(0 && "uncomputed matrix or unhandled weight type");
            return 0.0; /* warning suppressor */
            break;
    }

    return 0.0; /* warning suppressor */
}
double compute_distmatrix(instance inst) {
    if (inst->adjmatrix != NULL) return 0.0;

    /* initialize wall clock time */
    struct timespec s, e;
    s.tv_sec = e.tv_sec = -1;
    stopwatch(&s, &e);

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
                    assert(0 && "unhandeled distance function");
                    break;
            }
        }
    }

    return stopwatch(&s, &e);
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

/* graphs utils */
int reachable(int* succ, int i, int j) {
    int ans = 0;
    int next = i;
    do {
        next = succ[next];
        if (next == -1) return -1;

        ans++;

        if (next == j) return ans;
    } while (succ[next] != i);

    return -1;
}
int visitable(int* link, int nnodes) {
    int visits = 0;

    int next = 0;
    do {
        next = link[next];
        visits++;
    } while (link[next] != 0);

    return visits == nnodes - 1;
}
void reverse_path(int* succ, int nnodes, int start, int end) {
    int* stack = (int*)malloc(nnodes * sizeof(int));
    int k = 0;

    stack[0] = start;
    k++;

    while (stack[k - 1] != end) {
        if (stack[k - 1] == -1) {
            printf("mad\n");
        }
        stack[k] = succ[stack[k - 1]];
        k++;
    }

    k--;
    while (k >= 1) {
        succ[stack[k]] = stack[k - 1];
        k--;
    }

    free(stack);
}

/* edges representation convertes */
int* edges_tosucc(edge* edges, int nnodes) {
    /* return nonzero if not a tour, 0 otherwise */

    int* succ = (int*)malloc(nnodes * sizeof(int));
    int tovisit = nnodes;

    int i = 0;
    while (tovisit--) {
        int tofind = edges[i].j;

        int j = (i + 1) % nnodes;
        while (j != i) {
            /* if I found it in the second position, swap */
            if (edges[j].j == tofind) swap(&edges[j].i, &edges[j].j);
            if (edges[j].i == tofind) break;

            j = (j + 1) % nnodes;
        }

        if (j == i) {
            /* not found => tour not complete! */
            return NULL;
        }

        /* store the node in the succ and update i as j, so we can restart from
         * there: its a double not indexed cycle, a bit odd but functioning */
        succ[tofind] = edges[j].j;
        i = j;
    }

    return succ;
}

int* randomtour(int nnodes, unsigned int seedp) {
    int* succ = (int*)malloc(nnodes * sizeof(int));

    int left_nodes = nnodes;
    int* pool = (int*)malloc(nnodes * sizeof(int));
    for (int i = 0; i < nnodes; i++) pool[i] = i;

    int start, act, next;
    int pool_pick;

    pool_pick = rand_r(&seedp) % left_nodes;
    start = act = pool[pool_pick];
    swap(&pool[pool_pick], &pool[left_nodes - 1]);
    left_nodes--;

    while (left_nodes > 0) {
        pool_pick = rand_r(&seedp) % left_nodes;
        next = pool[pool_pick];
        swap(&pool[pool_pick], &pool[left_nodes - 1]);

        succ[act] = next;
        act = next;

        left_nodes--;
    }
    succ[next] = start;

    free(pool);

    return succ;
}

/* edge comparator for kruskal */
int wedgecmp(const void* a, const void* b) {
    double wa = ((wedge*)a)->w;
    double wb = ((wedge*)b)->w;

    return wa - wb < EPSILON ? -1 : +1;
}
/* compare nodes lexycographically */
int nodelexcmp(const void* a, const void* b) {
    node* na = (node*)a;
    node* nb = (node*)b;

    // TODO(lugot): CHECK a < b vs a - b < EPS FP SAFETY
    if (fabs(na->x - nb->x) > EPSILON) return na->x < nb->x ? -1 : 1;
    return na->y < nb->y ? -1 : 1;
}
int pathcmp(const void* a, const void* b, void* data) {
    int* pathlenghts = (int*)data;

    if (pathlenghts[*((int*)a)] < pathlenghts[*((int*)b)]) return -1;
    return 1;
}

/* computational geometry helpers */
double cross(node a, node b) { return a.x * b.y - a.y * b.x; }
int ccw(node a, node b, node c) {
    node ba = (node){b.x - a.x, b.y - a.y};
    node ca = (node){c.x - a.x, c.y - a.y};

    return cross(ba, ca) > EPSILON;
}

/* quick string helpers */
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
        case MTZ_LAZY_NO2:
            snprintf(ans, bufsize, "mtz_lazy_no2");
            break;
        case MTZ_INDICATOR:
            snprintf(ans, bufsize, "mtz_indicator");
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
        case SOFT_FIXING:
            snprintf(ans, bufsize, "soft_fixing");
            break;
        case MST:
            snprintf(ans, bufsize, "minspantree");
            break;
        case GREEDY:
            snprintf(ans, bufsize, "greedy");
            break;
        case GRASP:
            snprintf(ans, bufsize, "grasp");
            break;
        case EXTRA_MILEAGE:
            snprintf(ans, bufsize, "extra_mileage");
            break;
        case VNS:
            snprintf(ans, bufsize, "vns_randomstart");
            break;
        case TABU_SEACH:
            snprintf(ans, bufsize, "tabu_search_randomstart");
            break;
    }

    return ans;
}
char* model_folder_tostring(enum model_folders folder) {
    int bufsize = 100;
    char* ans = (char*)calloc(bufsize, sizeof(char));

    switch (folder) {
        case TSPLIB:
            snprintf(ans, bufsize, "tsplib");
            break;
        case GENERATED:
            snprintf(ans, bufsize, "generated");
            break;
    }

    return ans;
}

/* wall clock trackers */
int64_t stopwatch(struct timespec* s, struct timespec* e) {
    /* if stopwatch not started yet, do if */
    if (s->tv_sec == -1) clock_gettime(CLOCK_REALTIME, s);

    /* the next times do that for end and return difference
     * first call result can be ignored */
    clock_gettime(CLOCK_REALTIME, e);

    return (e->tv_sec - s->tv_sec) * 1000 +
           (time_t)((e->tv_nsec - s->tv_nsec)) / 1000000;
}
int64_t stopwatch_n(struct timespec* s, struct timespec* e) {
    /* if stopwatch not started yet, do if */
    if (s->tv_sec == -1) clock_gettime(CLOCK_REALTIME, s);

    /* the next times do that for end and return difference
     * first call result can be ignored */
    clock_gettime(CLOCK_REALTIME, e);

    return (e->tv_sec - s->tv_sec) * 1e9 + (time_t)(e->tv_nsec - s->tv_nsec);
}

/* classic utils */
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
void* intset(int* arr, int c, int n) {
    for (int i = 0; i < n; i++) arr[i] = c;

    return arr;
}
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
char** list_files(enum model_folders folder, int* nmodels) {
    char* path;
    switch (folder) {
        case TSPLIB:
            path = "../data_tsplib";
            break;
        case GENERATED:
            path = "../data_generated";
            break;
    }

    /* first count nfiles */
    DIR* dp = opendir(path);
    assert(dp != NULL);

    *nmodels = 0;
    struct dirent* ep;
    while ((ep = readdir(dp))) (*nmodels)++;
    closedir(dp);
    *(nmodels) -= 2;

    char** model_names = (char**)calloc(*nmodels, sizeof(char*));

    /* then grab em */
    dp = opendir(path);
    assert(dp != NULL);

    int i = 0;
    while ((ep = readdir(dp))) {
        if (!strcmp(ep->d_name, ".")) continue;
        if (!strcmp(ep->d_name, "..")) continue;

        model_names[i] = (char*)calloc(1 + strlen(ep->d_name), sizeof(char));
        snprintf(model_names[i], 1 + strlen(ep->d_name), "%s", ep->d_name);

        i++;
    }
    closedir(dp);

    return model_names;
}
