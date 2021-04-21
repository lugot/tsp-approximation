#include "../include/tsp.h"

#include <assert.h>
#include <time.h>

#include "../include/string.h"
#include "../include/utils.h"

cplex_params create_params() {
    cplex_params params =
        (cplex_params)calloc(1, sizeof(struct cplex_params_t));

    /* set defaults for params */
    params->num_threads = -1;
    params->timelimit = CPX_INFBOUND;
    params->available_memory = 4096;
    params->cost = REAL;

    return params;
}

instance create_empty_instance() {
    instance inst = (instance)calloc(1, sizeof(struct instance_t));
    inst->params = (cplex_params)calloc(1, sizeof(struct cplex_params_t));

    /* initializing only the non-zero default parameters */
    inst->params->num_threads = -1;
    inst->params->timelimit = CPX_INFBOUND;
    inst->params->available_memory = 4096;
    inst->params->cost = REAL;

    return inst;
}

instance create_instance(cplex_params params) {
    instance inst = (instance)calloc(1, sizeof(struct instance_t));

    /* initializing with passed parameters */
    inst->params = params;

    return inst;
}

void add_params(instance inst, cplex_params params) {
    /* should be non NULL by default constructor */
    memcpy(inst->params, params, sizeof(struct cplex_params_t));
}

instance generate_random_instance(int id, int nnodes, int box_size) {
    instance inst = create_empty_instance();
    srand(id);

    char* buf;
    int bufsize = 100;
    buf = (char*)calloc(bufsize, sizeof(char));

    snprintf(buf, bufsize, "random%d", id);
    inst->model_name = (char*)calloc(1 + strlen(buf), sizeof(char));
    snprintf(inst->model_name, bufsize, "%s", buf);

    snprintf(buf, bufsize, "random generated with id %d", id);
    inst->model_comment = (char*)calloc(1 + strlen(buf), sizeof(char));
    snprintf(inst->model_comment, bufsize, "%s", buf);

    inst->weight_type = ATT;
    inst->nnodes = nnodes;

    inst->nodes = (node*)calloc(nnodes, sizeof(struct node_t));
    for (int i = 0; i < nnodes; i++) {
        inst->nodes[i].x = (double)rand() / ((double)RAND_MAX / box_size);
        inst->nodes[i].y = (double)rand() / ((double)RAND_MAX / box_size);

        /* TODO: investigate
         *uint64_t r53 = ((uint64_t)(rand()) << 21) ^ (rand() >> 2);
         *return (double)r53 / 9007199254740991.0; // 2^53 - 1
         */
    }

    free(buf);

    return inst;
}

instance* generate_random_instances(int ninstances, int nnodes, int box_size) {
    instance* insts = (instance*)calloc(ninstances, sizeof(struct instance_t));

    // TODO(lugot): do it better
    srand(time(NULL));
    unsigned int seedp = 0;

    for (int i = 0; i < ninstances; i++) {
        /* between 0 and RAND_MAX (2^31 or smht here) */
        insts[i] = generate_random_instance(rand_r(&seedp), nnodes, box_size);
    }

    return insts;
}

// TODO(lugot): test
instance clone_instance(instance inst) {
    instance newone = (instance)calloc(1, sizeof(struct instance_t));

    newone->model_name =
        (char*)calloc(1 + strlen(inst->model_name), sizeof(char));
    snprintf(newone->model_name, 1 + strlen(inst->model_name), "%s",
             inst->model_name);
    newone->model_comment =
        (char*)calloc(1 + strlen(inst->model_comment), sizeof(char));
    snprintf(newone->model_comment, 1 + strlen(inst->model_comment), "%s",
             inst->model_comment);

    newone->params = (cplex_params)calloc(1, sizeof(struct cplex_params_t));
    memcpy(newone->params, inst->params, sizeof(struct cplex_params_t));

    return newone;
}

void free_instance(instance inst) {
    int nnodes = inst->nnodes;
    free(inst->model_name);
    free(inst->model_comment);

    free(inst->nodes);
    for (int i = 0; i < nnodes; i++) free(inst->adjmatrix[i]);
    free(inst->adjmatrix);

    // TODO(lugot): free sols?
    /*for (int i=0; i<nnodes; i++) free(inst->adjmatrix[i]);*/
    free(inst->sols);
}

void save_instance(instance inst) {
    char* filename;
    int bufsize = 100;
    filename = (char*)calloc(bufsize, sizeof(char));
    snprintf(filename, bufsize, "../data_generated/%s/%s.tsp", inst->model_name,
             inst->model_name);

    free(filename);
}

void add_solution(instance inst, solution sol) {
    int nsols = inst->nsols;
    inst->sols =
        (solution*)realloc(inst->sols, (nsols + 1) * sizeof(struct solution_t));
    inst->sols[nsols] = sol;
    sol->inst = inst;

    inst->nsols++;
}

/* printers */
void print_instance(instance inst, int print_data) {
    int nnodes = inst->nnodes;
    int nsols = inst->nsols;

    printf("\ninfos:\nmodel: %s\n", inst->model_name);
    printf("comment: %s\n", inst->model_comment);
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

    print_cplex_params(inst->params);

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
    printf("- number of nodes: %d\n", nnodes);
    if (print_data) {
        printf("- weights:\n");
        if (inst->adjmatrix != NULL) {
            /* compute numof figures for spacing */
            int column_width = 0;
            for (int i = 0; i < nnodes; i++)
                for (int j = 0; j <= i; j++) {
                    column_width = max(column_width,
                                       1 + (int)log10(inst->adjmatrix[i][j]));
                }
            column_width = 2 + max(column_width, 1 + (int)log10(nnodes));

            /* figures.ab so max 1 decimal figures */
            char* buf = (char*)calloc(1 + column_width, sizeof(char));

            for (int i = 0; i < nnodes; i++) {
                snprintf(buf, 1 + column_width, "%d | ", i + 1);
                printf("%*s", 1 + column_width, buf);

                for (int j = 0; j <= i; j++) {
                    snprintf(buf, 1 + column_width, "%.1lf ",
                             inst->adjmatrix[i][j]);
                    printf("%*s", 1 + column_width, buf);
                }
                printf("\n");
            }
            free(buf);
        }
        printf("- nodes:\n");
        if (inst->nodes != NULL) {
            int column_width = 0;
            for (int i = 0; i < nnodes; i++) {
                column_width =
                    max(column_width, 1 + (int)log10(inst->nodes[i].x));
                column_width =
                    max(column_width, 1 + (int)log10(inst->nodes[i].y));
            }
            column_width = 2 + max(column_width, 1 + (int)log10(nnodes));

            /* figures.ab so max 1 decimal figures */
            char* buf = (char*)calloc(1 + column_width, sizeof(char));

            for (int i = 0; i < nnodes; i++) {
                snprintf(buf, 1 + column_width, "%d | ", i + 1);
                printf("%*s", column_width + 1, buf);

                snprintf(buf, 1 + column_width, "%.1lf ", inst->nodes[i].x);
                printf("%*s", column_width + 1, buf);

                snprintf(buf, 1 + column_width, "%.1lf ", inst->nodes[i].y);
                printf("%*s\n", column_width + 1, buf);
            }
            free(buf);
        }
    }

    printf("solutions:\n");
    printf("- num of solutions: %d\n", nsols);
    for (int i = 0; i < nsols; i++) {
        printf("solution %d:\n", i + 1);
        print_solution(inst->sols[i], print_data);
    }

    printf("--- ---\n\n");
}

void print_cplex_params(cplex_params params) {
    printf("cplex params:\n");
    if (params == NULL) return;
    printf("- random seed: %d\n", params->randomseed);
    printf("- number of threads: %d\n", params->num_threads);
    printf("- time limit: %lf\n", params->timelimit);
    printf("- available memory: %d MB\n", params->available_memory);
    printf("- costs type: ");
    switch (params->cost) {
        case REAL:
            printf("real\n");
            break;
        case INTEGER:
            printf("integer\n");
            break;
    }
}

void print_solution(solution sol, int print_data) {
    int nedges = sol->nedges;

    printf("- optimality: ");
    switch (sol->model_type) {
        case OPTIMAL_TOUR:
            printf("optimal tour from file\n");
            break;
        case NOSEC:
            printf("symmetric, no subtour elimination\n");
            break;
        case MTZ_STATIC:
            printf("asymmetric, mtz subtour elimination\n");
            break;
        case MTZ_LAZY:
            printf("asymmetric, mtz subtour elimination, lazy constraints\n");
            break;
        case GGLIT_STATIC:
            printf("asymmetric, gg literature subtour elimination\n");
            break;
        case GGFISH_STATIC:
            printf("asymmetric, gg prof formulation subtour elimination\n");
            break;
        case GG_LAZY:
            printf("asymmetric, gg subtour elimination, lazy constraints\n");
            break;
        case BENDERS:
            printf("symmetric, benders loop method\n");
            break;
        case BENDERS_CALLBACK:
            printf("symmetric, benders loop method w/ callback\n");
            break;
    }
    printf("- zstar: %lf\n", sol->zstar);
    printf("- num edges: %d\n", nedges);
    if (print_data) {
        printf("- edges:\n");
        if (sol->edges != NULL) {
            int column_width = 2 + (int)log10(nedges);
            char* buf = (char*)calloc(column_width, sizeof(char));

            for (int i = 0; i < nedges; i++) {
                snprintf(buf, column_width, "%d | ", i + 1);
                printf("%*s", column_width + 2, buf);

                snprintf(buf, column_width, "%d ", sol->edges[i].i + 1);
                printf("%*s", column_width, buf);

                snprintf(buf, column_width, "%d ", sol->edges[i].j + 1);
                printf("%*s\n", column_width, buf);
            }
            free(buf);
        }
        printf("- parent:\n");
        if (sol->link != NULL) {
            int column_width = 2 + (int)log10(nedges);
            char* buf = (char*)calloc(column_width, sizeof(char));

            for (int i = 0; i < nedges; i++) {
                snprintf(buf, column_width, "%d", i + 1);
                printf("%*s", column_width + 1, buf);
            }
            printf("\n");
            for (int i = 0; i < nedges; i++) {
                snprintf(buf, column_width, "%d", sol->link[i] + 1);
                printf("%*s", column_width + 1, buf);
            }
            printf("\n");
            free(buf);
        }
    }
    printf("- distance time: %lf ms\n", sol->distance_time);
    printf("- build time: %lf ms\n", sol->build_time);
    printf("- solve time: %lf s\n", sol->solve_time / 1000.0);
}

void plot_solution_graphviz(solution sol) {
    plot_solutions_graphviz((solution[]){sol}, 1);
}

void plot_solutions_graphviz(solution* sols, int num_sols) {
    // TODO(lugot): make proportional to the number of nodes
    double box_size = 20.0;
    double max_coord = 0.0;

    instance inst = sols[0]->inst;
    assert(inst != NULL && "no instance associated with solution");
    int nnodes = inst->nnodes;

    for (int i = 0; i < nnodes && inst->nodes != NULL; i++) {
        max_coord = max(max_coord, fabs(inst->nodes[i].x));
        max_coord = max(max_coord, fabs(inst->nodes[i].y));
    }

    char* fname;
    int bufsize = 100;
    fname = (char*)calloc(bufsize, sizeof(char));
    if (inst->model_folder == TSPLIB)
        snprintf(fname, bufsize, "../data_tsplib/%s/%s.dot", inst->model_name,
                 inst->model_name);
    if (inst->model_folder == GENERATED)
        snprintf(fname, bufsize, "../data_generated/%s/%s.dot",
                 inst->model_name, inst->model_name);

    FILE* fp;
    fp = fopen(fname, "w");
    assert(fp != NULL && "file not found while saving .dot");

    fprintf(fp, "graph %s {\n", inst->model_name);
    fprintf(fp, "\tnode [shape=circle fillcolor=white]\n");
    for (int i = 0; i < nnodes && inst->nodes != NULL; i++) {
        double plot_posx = inst->nodes[i].x / max_coord * box_size;
        double plot_posy = inst->nodes[i].y / max_coord * box_size;

        fprintf(fp, "\t%d [ pos = \"%lf,%lf!\"]\n", i + 1, plot_posx,
                plot_posy);
    }
    fprintf(fp, "\n");

    for (int u = 0; u < num_sols; u++) {
        for (int k = 0; k < sols[u]->nedges; k++) {
            fprintf(fp, "\t%d -- %d", sols[u]->edges[k].i + 1,
                    sols[u]->edges[k].j + 1);

            if (sols[u]->model_type == OPTIMAL_TOUR)
                fprintf(fp, " [color = red]");
            fprintf(fp, "\n");
        }
    }
    fprintf(fp, "}");

    fclose(fp);
    free(fname);
}
