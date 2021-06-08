#include "../include/tsp.h"

#include <assert.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include "../include/string.h"
#include "../include/utils.h"
#include "globals.h"

/* cplex param manipulators */
cplex_params create_params() {
    cplex_params params = (cplex_params)malloc(sizeof(struct cplex_params_t));

    /* set defaults for params */
    params->randomseed = 1337;
    params->num_threads = -1;
    params->timelimit = CPX_INFBOUND;
    params->available_memory = 4096;
    params->cost = REAL;

    return params;
}
void add_params(instance inst, cplex_params params) {
    assert(inst != NULL);
    assert(params != NULL);

    inst->params->randomseed = params->randomseed;
    inst->params->num_threads = params->num_threads;
    inst->params->timelimit = params->timelimit;
    inst->params->available_memory = params->available_memory;
    inst->params->cost = params->cost;

    /* memcpy(inst->params, params, sizeof(struct cplex_params_t)); */
}
void free_params(cplex_params params) { free(params); }

/* instance manipulators */
instance create_empty_instance() {
    instance inst = (instance)calloc(1, sizeof(struct instance_t));
    inst->params = create_params();

    return inst;
}
instance create_instance(cplex_params params) {
    assert(params != NULL);

    instance inst = (instance)calloc(1, sizeof(struct instance_t));
    /* initializing with passed parameters */
    inst->params = params;

    return inst;
}
instance generate_random_instance(int id, int nnodes) {
    instance inst = create_empty_instance();
    srand(id);

    // random pick from an uniforme distribution the number of nodes in range
    // [9/10(nnodes), 11/10(nnodes)]
    double random_fluctuation = (double)rand() / (double)RAND_MAX - 0.5;
    random_fluctuation = random_fluctuation * nnodes / 5;
    nnodes += (int)random_fluctuation;

    printf("Number of nodes %d\n", nnodes);
    char* buf;
    int bufsize = 100;
    buf = (char*)calloc(bufsize, sizeof(char));

    snprintf(buf, bufsize, "random%d", id);
    inst->model_name = (char*)calloc(1 + strlen(buf), sizeof(char));
    snprintf(inst->model_name, bufsize, "%s", buf);

    snprintf(buf, bufsize, "random generated with id %d", id);
    inst->model_comment = (char*)calloc(1 + strlen(buf), sizeof(char));
    snprintf(inst->model_comment, bufsize, "%s", buf);

    inst->model_folder = GENERATED;
    inst->instance_type = TSP;

    inst->weight_type = ATT;
    inst->nnodes = nnodes;

    int box_size = 20.0;
    inst->nodes = (node*)calloc(nnodes, sizeof(struct node_t));
    for (int i = 0; i < nnodes; i++) {
        inst->nodes[i].x = (double)rand() / ((double)RAND_MAX / box_size);
        inst->nodes[i].y = (double)rand() / ((double)RAND_MAX / box_size);
    }

    free(buf);

    return inst;
}
instance* generate_random_instances(int ninstances, int nnodes) {
    instance* insts = (instance*)calloc(ninstances, sizeof(struct instance_t));

    // TODO(lugot): do it better
    srand(time(NULL));
    unsigned int seedp = 0;

    for (int i = 0; i < ninstances; i++) {
        /* between 0 and RAND_MAX (2^31 or smht here) */
        insts[i] = generate_random_instance(rand_r(&seedp), nnodes);
    }

    return insts;
}
void save_instance(instance inst) {
    assert(inst != NULL);

    int bufsize = 100;
    char* folder = model_folder_tostring(inst->model_folder);
    char* dirname = (char*)calloc(bufsize, sizeof(char));
    snprintf(dirname, bufsize, "../data_%s/%s", folder, inst->model_name);

    /* create the directory if not exists */
    struct stat st = {0};
    if (stat(dirname, &st) == -1) {
        mkdir(dirname, 0700);
    }

    /* create the file */
    char* fname;
    fname = (char*)calloc(bufsize, sizeof(char));

    snprintf(fname, bufsize, "../data_%s/%s/%s.tsp", folder, inst->model_name,
             inst->model_name);

    FILE* fp;
    fp = fopen(fname, "w");
    assert(fp != NULL && "file not found while saving .tsp");

    fprintf(fp, "NAME : %s\n", inst->model_name);
    fprintf(fp, "COMMENT: %s\n", inst->model_comment);
    fprintf(fp, "TYPE : TSP\n");
    fprintf(fp, "DIMENSION : %d\n", inst->nnodes);
    fprintf(fp, "EDGE_WEIGHT_TYPE : EUC_2D\n");
    fprintf(fp, "NODE_COORD_SECTION\n");
    for (int i = 0; i < inst->nnodes; i++) {
        fprintf(fp, "%d %f %f\n", i + 1, inst->nodes[i].x, inst->nodes[i].y);
    }
    fprintf(fp, "EOF\n");

    fclose(fp);
    free(fname);
    free(folder);
    free(dirname);
}
void free_instance(instance inst) {
    free(inst->model_name);
    free(inst->model_comment);

    free(inst->params);

    free(inst->nodes);
    if (inst->adjmatrix != NULL) {
        for (int i = 0; i < inst->nnodes; i++) free(inst->adjmatrix[i]);
        free(inst->adjmatrix);
    }

    for (int i = 0; i < inst->nsols; i++) free_solution(inst->sols[i]);
    free(inst->sols);

    free(inst);
}

/* solution manipulators */
solution create_solution(instance inst, enum model_types model_type,
                         int nedges) {
    solution sol = (solution)calloc(1, sizeof(struct solution_t));

    sol->model_type = model_type;
    sol->nedges = nedges;
    sol->edges = (edge*)calloc(nedges, sizeof(struct edge_t));

    return sol;
}
void add_solution(instance inst, solution sol) {
    int nsols = inst->nsols;
    inst->sols =
        (solution*)realloc(inst->sols, (nsols + 1) * sizeof(struct solution_t));
    inst->sols[nsols] = sol;
    sol->inst = inst;

    inst->nsols++;
}
void free_solution(solution sol) {
    char* model_type_str = model_type_tostring(sol->model_type);
    free(model_type_str);
    free(sol->edges);

    free(sol);
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
        if (inst->adjmatrix != NULL && EXTRA) {
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

    char* type = model_type_tostring(sol->model_type);
    printf("- model_type: %s\n", type);
    free(type);
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
    }
    printf("- distance time: %lf ms\n", sol->distance_time);
    printf("- build time: %lf ms\n", sol->build_time);
    printf("- solve time: %lf s\n", sol->solve_time / 1000.0);
}

/* plotters */
void plot_graphviz(solution sol, int* edgecolors, int version) {
    double box_size = 20.0;
    double max_coord = 0.0;

    instance inst = sol->inst;
    assert(inst != NULL && "no instance associated with solution");
    assert(inst->nodes != NULL);
    int nnodes = inst->nnodes;

    for (int i = 0; i < nnodes && inst->nodes != NULL; i++) {
        max_coord = max(max_coord, fabs(inst->nodes[i].x));
        max_coord = max(max_coord, fabs(inst->nodes[i].y));
    }

    int bsize = 100;
    char* fname = (char*)calloc(bsize, sizeof(char));
    char* folder = model_folder_tostring(inst->model_folder);

    snprintf(fname, bsize, "../data_%s/%s/%s.%d.dot", folder, inst->model_name,
             inst->model_name, version);

    free(folder);

    FILE* fp;
    fp = fopen(fname, "w");
    assert(fp != NULL && "file not found while saving .dot");

    fprintf(fp, "graph %s {\n", inst->model_name);
    if (inst->nnodes < 100) {
        fprintf(fp, "\tnode [shape=circle fillcolor=white]\n");
    } else {
        fprintf(fp, "\tnode [shape=point fillcolor=white]\n");
    }
    for (int i = 0; i < nnodes && inst->nodes != NULL; i++) {
        double plot_posx = inst->nodes[i].x / max_coord * box_size;
        double plot_posy = inst->nodes[i].y / max_coord * box_size;

        fprintf(fp, "\t%d [ pos = \"%lf,%lf!\"]\n", i + 1, plot_posx,
                plot_posy);
    }
    fprintf(fp, "\n");

    int ncolors = 5;
    char* colors[] = {"black", "red", "green", "blue", "purple"};

    for (int k = 0; k < sol->nedges; k++) {
        fprintf(fp, "\t%d -- %d", sol->edges[k].i + 1, sol->edges[k].j + 1);

        if (sol->model_type == OPTIMAL_TOUR) {
            fprintf(fp, " [color = red]");
        } else if (edgecolors != NULL) {
            fprintf(fp, " [color = %s]", colors[edgecolors[k] % ncolors]);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "}");

    /* TODO(lugot): FIX */
    /* bsize = 200; */
    /* char* command = (char*)calloc(bsize, sizeof(char)); */
    /* snprintf(command, bsize, */
    /*          "(cd ../data_%s/%s && dot -Kneato %s.%d.dot -Tpng > %s.%d.png)",
     */
    /*          model_folder_tostring(inst->model_folder), inst->model_name, */
    /*          inst->model_name, version, inst->model_name, version); */

    /* printf("command: |%s|\n", command); */
    /* system(command); */
    /* free(command); */

    fclose(fp);
    free(fname);
}
void plot_profiler(instance* insts, int ninstances) {
    assert(insts != NULL);
    assert(insts[0] != NULL);
    // ciao

    /* remove and create new fresh csv */
    remove("../results/results.csv");
    FILE* fp;
    fp = fopen("../results/results.csv", "w");
    assert(fp != NULL && "file not found while saving .csv");

    /* save the data */
    int nmodels = insts[0]->nsols;
    fprintf(fp, "%d,", nmodels);

    for (int i = 0; i < nmodels; i++) {
        enum model_types model_type = insts[0]->sols[i]->model_type;
        assert(model_type != NOSEC && model_type != OPTIMAL_TOUR);

        char* model_name_str = model_type_tostring(model_type);
        fprintf(fp, "%s", model_name_str);
        free(model_name_str);

        if (i < nmodels - 1)
            fprintf(fp, ",");
        else
            fprintf(fp, "\n");
    }

    for (int i = 0; i < ninstances; i++) {
        instance inst = insts[i];
        fprintf(fp, "%s,", inst->model_name);

        assert(inst->nsols == nmodels && "missing some solutions");

        for (int j = 0; j < nmodels; j++) {
            assert(inst->sols[j]->model_type == insts[0]->sols[j]->model_type);

            if (j < nmodels - 1)
                fprintf(fp, "%lf,", inst->sols[j]->solve_time);
            else
                fprintf(fp, "%lf\n", inst->sols[j]->solve_time);
        }
    }

    fclose(fp);

    /* generate the plot */
    // TODO(lugot): adjust timelimit
    system(
        "python3 ../results/perprof.py -D , -T 3600 -S 2 -M 2 "
        "../results/results.csv ../results/pp.pdf -P 'model comparisons'");
}
