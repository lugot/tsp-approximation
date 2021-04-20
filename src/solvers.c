#include "../include/solvers.h"

#include <assert.h>
#include <cplex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include "../include/globals.h"
#include "../include/model_builder.h"
#include "../include/tsp.h"
#include "../include/union_find.h"
#include "../include/utils.h"

solution TSPopt(instance inst, enum model_types model_type) {
    assert(inst != NULL);
    assert(inst->params != NULL && "no CPLEX params found");
    assert(inst->instance_type == TSP && "need TSP instance");

    /* create and populate solution */
    solution sol = (solution)calloc(1, sizeof(struct solution_t));
    sol->model_type = model_type;

    int nedges;
    nedges = sol->nedges = inst->nnodes;
    sol->edges = (edge*)calloc(nedges, sizeof(struct edge_t));
    sol->link = (int*)calloc(nedges, sizeof(int));

    sol->distance_time = compute_dist(inst);

    /* open CPLEX model */
    int error;
    CPXENVptr env = CPXopenCPLEX(&error);
    CPXLPptr lp = CPXcreateprob(env, &error, "TSP");

    /* set params:
     * - timelimit
     * - EpInteger: CPLEX tollerance to declare a variable integer, important in
     * bigM
     * - EpRightHandSide: the less or equal satisfied up to this tollerance,
     * usually 1e-5, with bigM 1e-9 */
    CPXsetdblparam(env, CPX_PARAM_TILIM, inst->params->timelimit);
    CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0);
    CPXsetdblparam(env, CPX_PARAM_EPRHS, 1e-9);

    /* populate enviorment with model data */
    sol->build_time = build_tsp_model(env, lp, inst, model_type);

    /* add callback if required */
    if (model_type == BENDERS) {
        CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE;
        if (CPXcallbacksetfunc(env, lp, contextid, add_BENDERS_sec_callback,
                               sol)) {
            print_error("CPXcallbacksetfunc() error");
        }
    }

    struct timeval start, end;
    gettimeofday(&start, NULL);

    /* solve! */
    if (CPXmipopt(env, lp)) print_error("CPXmipopt() error");

    /* store the optimal solution found by CPLEX */
    int ncols = CPXgetnumcols(env, lp);
    double* xstar = (double*)calloc(ncols, sizeof(double));
    if (CPXgetx(env, lp, xstar, 0, ncols - 1)) print_error("CPXgetx() error\n");

    switch (model_type) {
        case NOSEC:
        case BENDERS_CALLBACK:
            get_symmsol(xstar, nedges, sol->edges, sol->link);
            break;

        case BENDERS:
            get_symmsol(xstar, nedges, sol->edges, sol->link);

            while (!visitable(sol->link, nedges)) {
                /* add the constraints */
                add_BENDERS_sec(env, lp, sol);

                /* some time has passed, updat the timelimit */
                gettimeofday(&end, NULL);
                int restime =
                    inst->params->timelimit - (end.tv_sec - start.tv_sec);
                CPXsetdblparam(env, CPX_PARAM_TILIM, restime);

                /* solve! */
                if (CPXmipopt(env, lp)) print_error("CPXmipopt() error");

                /* store the optimal solution found by CPLEX */
                int num_cols = CPXgetnumcols(env, lp);

                free(xstar);
                xstar = (double*)calloc(num_cols, sizeof(double));
                if (CPXgetx(env, lp, xstar, 0, num_cols - 1))
                    print_error("CPXgetx() error\n");

                /* retrive the solution */
                get_symmsol(xstar, nedges, sol->edges, sol->link);
                free(xstar);
            }

            /* save the complete model */
            char* filename;
            filename = (char*)calloc(100, sizeof(char));
            if (inst->model_folder == TSPLIB)
                snprintf(filename, 28 + 2 * strlen(inst->model_name),
                         "../data_tsplib/%s/%s.benders.lp", inst->model_name,
                         inst->model_name);
            if (inst->model_folder == GENERATED)
                snprintf(filename, 31 + 2 * strlen(inst->model_name),
                         "../data_generated/%s/%s.benders.lp", inst->model_name,
                         inst->model_name);
            CPXwriteprob(env, lp, filename, NULL);
            free(filename);

            break;

        case MTZ_STATIC:
        case MTZ_LAZY:
        case GG_STATIC:
        case GG_LAZY:
            get_asymmsol(xstar, nedges, sol->edges, sol->link);
            break;

        case OPTIMAL_TOUR:
            assert(model_type != OPTIMAL_TOUR &&
                   "tried to solve an optimal tour instance");
            break;
    }

    free(xstar);

    gettimeofday(&end, NULL);
    sol->solve_time = end.tv_sec - start.tv_sec;

    /* retreive the min cost */
    CPXgetobjval(env, lp, &sol->zstar);
    /*sol->zstar = zstar(inst, sol);*/

    /* add the solution to the pool associated with it's instance */
    add_solution(inst, sol);

    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);

    return sol;
}

void get_symmsol(double* xstar, int nedges, edge* edges, int* link) {
    assert(link != NULL);

    /* set up parent as linked list, no check on size */
    for (int i = 0; i < nedges; i++) link[i] = i;

    /* index over selected edges */
    int k = 0;

    /* check for edges j>i if is selected */
    for (int i = 0; i < nedges; i++) {
        for (int j = i + 1; j < nedges; j++) {
            if (xstar[xpos(i, j, nedges)] > 0.5) {
                /* link is mandatory, edges is optional */
                if (edges != NULL) edges[k++] = (edge){i, j};

                /* actually a shifted linked list */
                if (!reachable(link, i, j) && !reachable(link, j, i))
                    swap(&link[i], &link[j]);
            }
        }
    }

    assert(k == nedges && "not enought edges CPLEX solution");
}

void get_asymmsol(double* xstar, int nedges, edge* edges, int* link) {
    assert(link != NULL);

    /* set up parent as linked list */
    for (int i = 0; i < nedges; i++) link[i] = i;

    /* index over selected edges */
    int k = 0;

    /* check for all edges if is selected */
    for (int i = 0; i < nedges; i++) {
        for (int j = 0; j < nedges; j++) {
            if (xstar[xxpos(i, j, nedges)] > 0.5) {
                /* link is mandatory, edges is optional */
                if (edges != NULL) edges[k++] = (edge){i, j};

                /* actually a shifted linked list */
                if (!reachable(link, i, j)) swap(&link[i], &link[j]);
            }
        }
    }

    assert(k == nedges && "not enought edges CPLEX solution");
}

void save_results(instance* insts, int ninstances) {
    assert(insts != NULL);
    assert(insts[0] != NULL);

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
        "python ../results/perprof.py -D , -T 3600 -S 2 -M 20 "
        "../results/results.csv ../results/pp.pdf -P 'model comparisons'");
}
