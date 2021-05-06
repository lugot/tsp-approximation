#include "../include/solvers.h"

#include <assert.h>
#include <cplex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include "../include/adjlist.h"
#include "../include/globals.h"
#include "../include/model_builder.h"
#include "../include/tsp.h"
#include "../include/union_find.h"
#include "../include/utils.h"

void assert_correctness(solution sol);
void perform_BENDERS(CPXENVptr env, CPXLPptr lp, instance inst, double* xstar,
                     struct timespec s, struct timespec e);
void perform_HARD_FIXING(CPXENVptr env, CPXLPptr lp, instance inst,
                         double* xstar, int perc);
void perform_SOFT_FIXING(CPXENVptr env, CPXLPptr lp, instance inst,
                         double* xstar);

solution TSPopt(instance inst, enum model_types model_type) {
    assert(inst != NULL);
    assert(inst->params != NULL && "no CPLEX params found");
    assert(inst->instance_type == TSP && "need TSP instance");

    /* create and populate solution */
    solution sol = create_solution(inst, model_type, inst->nnodes);
    sol->distance_time = compute_dist(inst);
    /* sol->timetype = inst->timetype; */

    /* open CPLEX model */
    int error;
    CPXENVptr env = CPXopenCPLEX(&error);
    CPXLPptr lp = CPXcreateprob(env, &error, "TSP");

    CPXgetdettime(env, &sol->start);

    /* set params:
     * - timelimit (timetype = 0 -> seconds, timetype = 1 -> ticks)
     * - EpInteger: CPLEX tollerance to declare a variable integer, important in
     * bigM
     * - EpRightHandSide: the less or equal satisfied up to this tollerance,
     * usually 1e-5, with bigM 1e-9 */
    CPXsetdblparam(env, CPX_PARAM_TILIM, inst->params->timelimit);
    CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0);
    CPXsetdblparam(env, CPX_PARAM_EPRHS, 1e-9);
    // CPXsetintparam(env, CPX_PARAM_RANDOMSEED, seed); todo: add a seed to
    // Cplexrun (Passed by TSP opt?) TODO

    /* save the log of the execution */
    char logfile[] = "execution.log";
    CPXsetlogfilename(env, logfile, "w");

    /* populate enviorment with model data */
    sol->build_time = build_tsp_model(env, lp, inst, model_type);

    /* model preprocessing: add some callbacks or set params before execution */
    switch (model_type) {
        case BENDERS_CALLBACK: {
            CPXLONG contextid =
                CPX_CALLBACKCONTEXT_CANDIDATE | CPX_CALLBACKCONTEXT_RELAXATION;
            if (CPXcallbacksetfunc(env, lp, contextid,
                                   add_BENDERS_sec_callback_driver, sol)) {
                print_error("CPXcallbacksetfunc() error");
            }
            break;
        }

        case HARD_FIXING: {
            /* set iniitial fraction of timelimit for first execution */
            CPXsetdblparam(env, CPX_PARAM_TILIM,
                           inst->params->timelimit / HF_INITIAL_FRACTION_TIME);
            /* alternative: work on branching node limit */
            /* CPXsetintparam(env, CPX_PARAM_NODELIM, 10); */

            /* recompute timelimit considering time spent on intial solution */
            inst->params->timelimit *= (1 - 1.0 / HF_INITIAL_FRACTION_TIME);

            /* set benders callbacks as blackbox for matheuristic */
            CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE;
            if (CPXcallbacksetfunc(env, lp, contextid,
                                   add_BENDERS_sec_callback_driver, sol)) {
                print_error("CPXcallbacksetfunc() error");
            }
            break;
        }

        case SOFT_FIXING: {
            /* set iniitial fraction of timelimit for first execution */
            /* CPXsetdblparam(env, CPX_PARAM_TILIM, */
            /*                inst->params->timelimit /
             * SF_INITIAL_FRACTION_TIME); */
            /* alternative: work on branching node limit */
            CPXsetintparam(env, CPX_PARAM_NODELIM, 10);

            /* recompute timelimit considering time spent on intial solution */
            inst->params->timelimit *= (1 - 1.0 / SF_INITIAL_FRACTION_TIME);

            /* set benders callbacks as blackbox for matheuristic */
            CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE;
            if (CPXcallbacksetfunc(env, lp, contextid,
                                   add_BENDERS_sec_callback_driver, sol)) {
                print_error("CPXcallbacksetfunc() error");
            }
            break;
        }

        default:
            /* no preprocessing steps */
            break;
    }

    /* initialize total wall-clock time of execution */
    struct timespec s, e;
    s.tv_sec = e.tv_sec = -1;
    stopwatch(&s, &e);

    /* solve! */
    if (CPXmipopt(env, lp)) print_error("CPXmipopt() error");

    /* store the optimal solution found by CPLEX */
    int ncols = CPXgetnumcols(env, lp);
    double* xstar = (double*)calloc(ncols, sizeof(double));
    if (CPXgetx(env, lp, xstar, 0, ncols - 1)) print_error("CPXgetx() error\n");

    /* retreive solution or perform future steps */
    switch (model_type) {
        case NOSEC:
        case BENDERS_CALLBACK:
            get_symmsol(xstar, sol);
            break;

        case BENDERS:
            perform_BENDERS(env, lp, inst, xstar, s, e);
            get_symmsol(xstar, sol);
            break;

        case HARD_FIXING:
            perform_HARD_FIXING(env, lp, inst, xstar, HF_PERCENTAGE);
            get_symmsol(xstar, sol);
            break;

        case SOFT_FIXING:
            perform_SOFT_FIXING(env, lp, inst, xstar);
            get_symmsol(xstar, sol);
            break;

        case MTZ_STATIC:
        case MTZ_LAZY:
        case GGLIT_STATIC:
        case GGLECT_STATIC:
        case GGLIT_LAZY:
            get_asymmsol(xstar, sol);
            break;

        case OPTIMAL_TOUR:
            assert(model_type != OPTIMAL_TOUR &&
                   "tried to solve an optimal tour instance");
            break;
    }

    free(xstar);

    /* retreive the time, in ticks and wall-clock */
    CPXgetdettime(env, &sol->end);
    sol->solve_time = stopwatch(&s, &e);

    /* retreive the min cost */
    CPXgetobjval(env, lp, &sol->zstar);
    /*sol->zstar = compute_zstar(inst, sol);*/

    /* assert_correctness(sol); // TODO(any): waiting for rewriting */

    /* add the solution to the pool associated with it's instance */
    add_solution(inst, sol);

    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);

    return sol;
}

void perform_BENDERS(CPXENVptr env, CPXLPptr lp, instance inst, double* xstar,
                     struct timespec s, struct timespec e) {
    int nedges = inst->nnodes;

    /* create an adjacency list to track the edges .. */
    adjlist l = adjlist_create(nedges);
    /* .. and first fill it with the first solution found (with possible
     * subtours) */
    for (int i = 0; i < nedges; i++) {
        for (int j = i + 1; j < nedges; j++) {
            if (xstar[xpos(i, j, nedges)] > 0.5) {
                adjlist_add_edge(l, i, j);
            }
        }
    }

    /* iterate and add SEC until we obtain a singe tour */
    while (!adjlist_single_tour(l)) {
        /* add the constraints */
        add_BENDERS_sec(env, lp, l);

        /* some time has passed, update the timelimit the timelimit
         * we have to pass time in second! */
        long restime = inst->params->timelimit - stopwatch(&s, &e) / 1000.0;
        CPXsetdblparam(env, CPX_PARAM_TILIM, restime);

        /* solve again! */
        if (CPXmipopt(env, lp)) print_error("CPXmipopt() error");

        /* store the optimal solution found by CPLEX */
        if (CPXgetx(env, lp, xstar, 0, CPXgetnumcols(env, lp) - 1)) {
            print_error("CPXgetx() error\n");
        }

        /* completely reset the adjacency list without recreating the object and
         * fill it with the found solution */
        adjlist_hard_reset(l);
        for (int i = 0; i < nedges; i++) {
            for (int j = i + 1; j < nedges; j++) {
                if (xstar[xpos(i, j, nedges)] > 0.5) {
                    adjlist_add_edge(l, i, j);
                }
            }
        }
    }

    /* save the now complete model */
    char* filename;
    int bufsize = 100;
    filename = (char*)calloc(bufsize, sizeof(char));
    char* folder = model_folder_tostring(inst->model_folder);

    snprintf(filename, bufsize, "../data_%s/%s/%s.benders.lp", folder,
             inst->model_name, inst->model_name);

    CPXwriteprob(env, lp, filename, NULL);
    free(folder);
    free(filename);
}

void perform_HARD_FIXING(CPXENVptr env, CPXLPptr lp, instance inst,
                         double* xstar, int perc) {
    assert(perc >= 0 && perc < 100);

    int nedges = inst->nnodes;
    int ncols = inst->ncols;

    /* TODO(lugot): FIX */
    srand(0);
    unsigned int seedp = 0L;

    /* objective tracking */
    double best_obj, prev_obj;
    CPXgetobjval(env, lp, &prev_obj);

    /* preparing some constants to quickly fix bounds */
    const char lbc = 'L';
    const char ubc = 'U';
    const double one = 1.0;
    const double zero = 0.0;

    /* perform HARD FIXING:
     * for each iteration fix randomly some edges, solve the model and release
     * the nodes fixed in order to warm start the next iteration */
    for (int iter = 0; iter < HF_ITERATIONS; iter++) {
        /* create an adjacency to track the fixed edges */
        adjlist l = adjlist_create(nedges);

        /* iterate over edges to randomly fix them
         * not a single for cause of arc tracking */
        for (int i = 0; i < nedges; i++) {
            for (int j = i + 1; j < nedges; j++) {
                int pos = xpos(i, j, nedges);

                /* check if the edge is an actual edge and select it */
                if (xstar[pos] > 0.5 && rand_r(&seedp) % 100 < perc) {
                    /* if (VERBOSE) { */
                    /*     printf("[VERBOSE] fixing arc (%d,%d)\n", i + 1, j +
                     * 1); */
                    /* } */

                    /* fix the edge by moving the lower bound to one */
                    CPXchgbds(env, lp, 1, &pos, &lbc, &one);

                    /* track which edges has been selected */
                    adjlist_add_edge(l, i, j);
                }
            }
        }

        /* iterate over tracked edges to retreive the loose ends:
         * let's create a simple SEC by fix that edge to zero */
        int i, j;
        while (adjlist_get_loose_ends(l, &i, &j)) {
            int pos = xpos(i, j, nedges);

            /* check if the edge is not part of a path (alone): in that case we
             * have to keep the lower bound to one to maintain consistency */
            double lb;
            CPXgetlb(env, lp, &lb, pos, pos);
            /* if lb is 1.0 means that we previously fix it -> skip */
            if (lb == 1.0) continue;

            /* if (VERBOSE) printf("[VERBOSE] removing (%d,%d)\n", i + 1, j +
             * 1); */
            CPXchgbds(env, lp, 1, &pos, &ubc, &zero);
        }

        /* save the model for analysis */
        if (EXTRA) CPXwriteprob(env, lp, "./hard_fixing", "LP");

        /* reset the timelimit: fraction of number of iterations */
        CPXsetdblparam(env, CPX_PARAM_TILIM,
                       inst->params->timelimit / HF_ITERATIONS);
        /* alternative: set the nodelimit */
        /* CPXsetintparam(env, CPX_PARAM_NODELIM, 100); */

        /* solve! and get solution */
        if (CPXmipopt(env, lp)) print_error("CPXmipopt() error");
        if (CPXgetx(env, lp, xstar, 0, ncols - 1)) {
            print_error("CPXgetx() error\n");
        }

        if (VERBOSE) {
            /* store the objective */
            CPXgetobjval(env, lp, &best_obj);

            printf(
                "[VERBOSE]: hard fixing status (iter %d out of %d)\n"
                "\tprev_obj: %.20lf\n"
                "\t act_obj: %.20lf\n",
                iter, HF_ITERATIONS, prev_obj, best_obj);

            prev_obj = best_obj;
        }

        if (EXTRA) {
            /* create a solution for the plot */
            solution sol = create_solution(inst, HARD_FIXING, nedges);
            sol->inst = inst;
            get_symmsol(xstar, sol);

            /* track which edges we fixed */
            int* edgecolors = (int*)calloc(nedges, sizeof(int));
            int u, v;
            adjlist_reset(l);
            while (adjlist_get_edge(l, &u, &v)) {
                for (int z = 0; z < nedges; z++) {
                    if (sol->edges[z].i == z && sol->edges[z].j == v) {
                        edgecolors[z] = 1;
                        break;
                    }
                }
            }

            /* save the plot */
            /* plot_solution_graphviz(sol, edgecolors, iter); TODO(any): FIX for
             * generated */
        }

        /* free the adjlist, do it not for additional plot option */
        adjlist_free(l);

        /* relax the fixing: no need to check because most of the nodes are
         * fixed. Just checking if this is not the last iteration to provide a
         * coeherent model */
        if (iter != HF_ITERATIONS - 1) {
            for (int col = 0; col < ncols; col++) {
                CPXchgbds(env, lp, 1, &col, &lbc, &zero);
                CPXchgbds(env, lp, 1, &col, &ubc, &one);
            }
        }
    }

    if (VERBOSE) {
        printf(
            "[VERBOSE]: hard fixing status (on exit)\n"
            "\tprev_obj: %.20lf\n"
            "\t act_obj: %.20lf\n",
            prev_obj, best_obj);
    }
}

void perform_SOFT_FIXING(CPXENVptr env, CPXLPptr lp, instance inst,
                         double* xstar) {
    int nedges = inst->nnodes;
    int ncols = inst->ncols;

    double k = SF_INITIAL_K;
    double best_obj, prev_obj;
    CPXgetobjval(env, lp, &prev_obj);

    while (inst->params->timelimit > 0) {
        /* track the time for the current iteration */
        struct timespec s, e;
        s.tv_sec = e.tv_sec = -1;
        stopwatch(&s, &e);

        /* set the new timelimit and perform another resolution */
        CPXsetdblparam(env, CPX_PARAM_TILIM, inst->params->timelimit);

        char** cname = (char**)calloc(1, sizeof(char*));
        cname[0] = (char*)calloc(100, sizeof(char));

        // TODO(any): comment
        const char sense = 'L';
        const double rhs = k - nedges;

        /* fetch new row .. */
        int lastrow = CPXgetnumrows(env, lp);
        snprintf(cname[0], strlen(cname[0]), "soft_fixing");
        /* .. and fix the sense and rhs */
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
            print_error("wrong CPXnewrows [soft fixing]");
        }

        /* TODO(lugot): MODIFY double for for tracking fixed nodes */
        for (int i = 0; i < ncols; i++) {
            if (xstar[i] > 0.5) {
                /* change the coefficent from 0 to -1.0 for variables = 1 in
                 * xstar */
                if (CPXchgcoef(env, lp, lastrow, i, -1.0)) {
                    print_error("wrong CPXchgcoef [soft fixing]");
                }
            }
        }

        /* solve! and get solution */
        if (CPXmipopt(env, lp)) print_error("CPXmipopt() error");
        if (CPXgetx(env, lp, xstar, 0, ncols - 1)) {
            print_error("CPXgetx() error\n");
        }

        /* compute objective and compare it with the previous one to determine
         * if we need to enlarge the neighborhood size */
        CPXgetobjval(env, lp, &best_obj);

        if (VERBOSE) {
            printf(
                "[VERBOSE]: soft fixing status (bf k updating)\n"
                "\tremaining_time: %lf\n"
                "\tk: %lf\n"
                "\tprev_obj: %.20lf\n"
                "\t act_obj: %.20lf\n",
                inst->params->timelimit, k, prev_obj, best_obj);
        }

        /* if no improvment enlarge neigborhood size */
        if (fabs(prev_obj - best_obj) < EPSILON) k += SF_K_STEP;
        /* cut the computation if k too high */
        if (k >= SF_MAX_K) break;
        /* and update the objective */
        prev_obj = best_obj;

        /* update timelimit */
        inst->params->timelimit -= stopwatch(&s, &e) / 1000.0;

        if (inst->params->timelimit > 0)
            if (CPXdelrows(env, lp, lastrow, lastrow))
                print_error("CPXdelrows() error\n");
    }

    if (VERBOSE) {
        printf(
            "[VERBOSE]: soft fixing status (on exiting)\n"
            "\tremaining_time: %lf\n"
            "\tk: %lf\n"
            "\t act_obj: %.20lf\n",
            inst->params->timelimit, k, best_obj);
    }
}

void get_symmsol(double* xstar, solution sol) {
    /* index over selected edges */
    int k = 0;

    /* check for edges j>i if is selected */
    for (int i = 0; i < sol->nedges; i++) {
        for (int j = i + 1; j < sol->nedges; j++) {
            if (xstar[xpos(i, j, sol->nedges)] > 0.5) {
                /* then store it in the solution*/
                sol->edges[k] = (edge){i, j};

                k++;
            }
        }
    }

    // TODO(any): check this
    /*assert(k == nedges && "not enought edges CPLEX solution");*/
}

void get_asymmsol(double* xstar, solution sol) {
    /* index over selected edges */
    int k = 0;

    /* check for edges i, j if is selected */
    for (int i = 0; i < sol->nedges; i++) {
        for (int j = 0; j < sol->nedges; j++) {
            if (xstar[xpos(i, j, sol->nedges)] > 0.5) {
                /* then store it in the solution*/
                sol->edges[k] = (edge){i, j};

                k++;
            }
        }
    }

    // TODO(any): check this
    /*assert(k == nedges && "not enought edges CPLEX solution");*/
}

/* void assert_correctness(solution sol) { */
/*     TODO(any): rewrite why not */
/*     int visited = 0; */
/*     int act = 0; */
/*     do { */
/*         visited++; */
/*     } while ((act = sol->link[act]) != 0); */

/*     assert(visited == sol->nedges); */
/* } */
