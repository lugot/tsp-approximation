#include "../include/solvers.h"

#include <assert.h>
#include "../include/globals.h"
#include "../include/model_builder.h"
#include "../include/models/benders.h"
#include "../include/models/fixing.h"
#include "../include/utils.h"
#include "../include/approximations.h"
#include "../include/constructives.h"
#include "../include/refinements.h"
#include "../include/metaheuristics.h"

void assert_correctness(solution sol);
solution TSPopt(instance inst, enum model_types model_type);

solution solve(instance inst, enum model_types model_type) {
    /* initialize total wall-clock time of execution */
    struct timespec s, e;
    s.tv_sec = e.tv_sec = -1;
    stopwatch(&s, &e);

    solution sol;

    switch (model_type) {
        case NOSEC:
        case MTZ_STATIC:
        case MTZ_LAZY:
        case MTZ_LAZY_DEG2:
        case MTZ_LAZY_DEG3:
        case MTZ_INDICATOR:
        case GGLIT_STATIC:
        case GGLECT_STATIC:
        case GGLIT_LAZY:
        case GGLECT_LAZY:
        case GGLIT_STATIC_DEG2:
        case BENDERS:
        case BENDERS_TWOPHASES:
        case BENDERS_CALLBACK:
        case HARD_FIXING:
        case SOFT_FIXING:
            sol = TSPopt(inst, model_type);
            break;

        case MST:
            sol = TSPminspantree(inst);
            break;

        case GREEDY:
            sol = TSPgreedy(inst);
            break;

        case GRASP:
            sol = TSPgrasp(inst, 0);
            break;

        case EXTRA_MILEAGE:
            sol = TSPextramileage(inst);
            break;

        case TWOOPT_MULTISTART:
            sol = NULL;
            // TODO(lugot): IMPLEMENT
            break;

        case VNS_RANDOM:
            sol = TSPvns(inst, NULL);
            break;

        case VNS_GREEDY:
            sol = TSPvns(inst, NULL);
            break;

        case TABU_SEACH:
            sol = TSPtabusearch(inst, NULL);
            break;

        case OPTIMAL_TOUR:
            sol = NULL;
            assert(model_type != OPTIMAL_TOUR &&
                   "tried to solve an optimal tour instance");
            break;
    }
    sol->solve_time = stopwatch(&s, &e);

    return sol;
}

solution TSPopt(instance inst, enum model_types model_type) {
    assert(inst != NULL);
    assert(inst->params != NULL && "no CPLEX params found");

    /* create and populate solution */
    solution sol = create_solution(inst, model_type, inst->nnodes);

    /* open CPLEX model */
    int error;
    CPXENVptr env = CPXopenCPLEX(&error);
    CPXLPptr lp = CPXcreateprob(env, &error, "TSP");

    /* set params:
     * - timelimit (timetype = 0 -> seconds, timetype = 1 -> ticks)
     * - EpInteger: CPLEX tollerance to declare a variable integer, important in
     * bigM
     * - EpRightHandSide: the less or equal satisfied up to this tollerance,
     * usually 1e-5, with bigM 1e-9 */
    CPXsetdblparam(env, CPX_PARAM_TILIM, inst->params->timelimit);
    double epint = 0.0; /*  suppress warning on lazy constraints GG */
    if (model_type == GGLIT_LAZY || model_type == GGLECT_LAZY ||
        model_type == MTZ_LAZY || model_type == MTZ_LAZY_DEG2 ||
        model_type == MTZ_LAZY_DEG3) {
        /* double the EPRHS, but even with bigM, this is safe because our M is
         * O(n) */
        epint = 2e-9;
    }
    CPXsetdblparam(env, CPX_PARAM_EPINT, epint);
    CPXsetdblparam(env, CPX_PARAM_EPRHS, 1e-9);

    // CPXsetintparam(env, CPX_PARAM_RANDOMSEED, seed); todo: add a seed to
    // Cplexrun (Passed by TSP opt?) TODO

    /* save the log of the execution */
    int bufsize = 100;
    char* logfile = (char*)malloc(bufsize * sizeof(char));
    char* model_type_str = model_type_tostring(model_type);
    snprintf(logfile, bufsize, "execution_%s_%s.log", inst->instance_name,
             model_type_str);
    CPXsetlogfilename(env, logfile, "a");
    free(model_type_str);
    free(logfile);

    /* populate enviorment with model data */
    sol->build_time = build_tsp_model(env, lp, inst, model_type);

    /* model preprocessing: add some callbacks or set params before execution */
    switch (model_type) {
        case BENDERS_TWOPHASES: {
            CPXsetdblparam(env, CPX_PARAM_NODELIM, BENDERS2P_NODELIM);
        } break;

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
                           inst->params->timelimit * HF_INITIAL_PERC_TIME);
            /* alternative: work on branching node limit */
            /* CPXsetintparam(env, CPX_PARAM_NODELIM, 10); */

            /* recompute timelimit considering time spent on intial solution */
            inst->params->timelimit *= (1 - HF_INITIAL_PERC_TIME);

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
            CPXsetdblparam(env, CPX_PARAM_TILIM,
                           inst->params->timelimit * SF_INITIAL_PERC_TIME);
            /* alternative: work on branching node limit */
            CPXsetintparam(env, CPX_PARAM_NODELIM, 10);

            /* recompute timelimit considering time spent on intial solution */
            inst->params->timelimit *= (1 - SF_INITIAL_PERC_TIME);

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

    /* track deterministic time */
    CPXgetdettime(env, &sol->start);

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

        case BENDERS: {
            struct timespec s, e;
            s.tv_sec = e.tv_sec = -1;
            stopwatch(&s, &e);
            perform_BENDERS(env, lp, inst, xstar, s, e, 2);
            get_symmsol(xstar, sol);
            break;
        }

        case BENDERS_TWOPHASES: {
            struct timespec s, e;
            s.tv_sec = e.tv_sec = -1;
            stopwatch(&s, &e);
            perform_BENDERS(env, lp, inst, xstar, s, e, 1);
            get_symmsol(xstar, sol);
            break;
        }

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
        case MTZ_LAZY_DEG2:
        case MTZ_LAZY_DEG3:
        case MTZ_INDICATOR:
        case GGLIT_STATIC:
        case GGLECT_STATIC:
        case GGLIT_LAZY:
        case GGLECT_LAZY:
        case GGLIT_STATIC_DEG2:
            get_asymmsol(xstar, sol);
            break;

        case OPTIMAL_TOUR:
            assert(model_type != OPTIMAL_TOUR &&
                   "tried to solve an optimal tour instance");
            break;

        case MST:
        case GREEDY:
        case GRASP:
        case EXTRA_MILEAGE:
        case TWOOPT_MULTISTART:
        case VNS_RANDOM:
        case VNS_GREEDY:
        case TABU_SEACH:
            assert(0 == 1 && "tried to solve a metaheuristic");
            break;
    }

    free(xstar);

    /* retreive the time, in ticks and wall-clock */
    CPXgetdettime(env, &sol->end);

    /* retreive the min cost */
    CPXgetobjval(env, lp, &sol->zstar);
    /*sol->zstar = compute_zstar(inst, sol);*/

    /* char solfile[] = "solution.xml"; */
    /* CPXsolwrite(env, lp, solfile); */

    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);

    return sol;
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
            if (xstar[xxpos(i, j, sol->nedges)] > 0.5) {
                /* then store it in the solution*/
                sol->edges[k] = (edge){i, j};

                k++;
            }
        }
    }

    // TODO(any): check this
    /*assert(k == nedges && "not enought edges CPLEX solution");*/
}
