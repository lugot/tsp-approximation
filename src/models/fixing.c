#include "../../include/models/fixing.h"

#include <assert.h>
#include <string.h>

#include "../../include/adjlist.h"
#include "../../include/globals.h"
#include "../../include/solvers.h"
#include "../../include/utils.h"

void perform_HARD_FIXING(CPXENVptr env, CPXLPptr lp, instance inst,
                         double* xstar, int perc) {
    assert(perc >= 0 && perc < 100);

    int nedges = inst->nnodes;
    int ncols = inst->ncols;

    double original_timelimit = inst->params->timelimit;
    inst->params->timelimit *= (1 - HF_INITIAL_PERC_TIME);

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
    int iter = 0;
    while (inst->params->timelimit > 0) {
        struct timespec s, e;
        s.tv_sec = e.tv_sec = -1;
        stopwatch(&s, &e);

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
        if (EXTRA_VERBOSE) CPXwriteprob(env, lp, "./hard_fixing", "LP");

        /* reset the timelimit: fraction of number of iterations */
        CPXsetdblparam(env, CPX_PARAM_TILIM,
                       inst->params->timelimit);
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
                "[VERBOSE]: hard fixing status (timelimit %lf)\n"
                "\tprev_obj: %.20lf\n"
                "\t act_obj: %.20lf\n",
                inst->params->timelimit, prev_obj, best_obj);

            prev_obj = best_obj;
        }

        if (EXTRA_VERBOSE) {
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
            plot_graphviz(sol, edgecolors, iter);
        }

        /* free the adjlist, do it not for additional plot option */
        adjlist_free(l);

        /* update timelimit */
        inst->params->timelimit -= stopwatch(&s, &e) / 1000.0;
        iter++;
        /* relax the fixing: no need to check because most of the nodes are
         * fixed. Just checking if this is not the last iteration to provide a
         * coeherent model */
        if (inst->params->timelimit > 0) {
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

    inst->params->timelimit = original_timelimit;
}

void perform_SOFT_FIXING(CPXENVptr env, CPXLPptr lp, instance inst,
                         double* xstar) {
    int nedges = inst->nnodes;
    int ncols = inst->ncols;


    double original_timelimit = inst->params->timelimit;
    inst->params->timelimit *= (1 - SF_INITIAL_PERC_TIME);
    double iter_time_limit = inst->params->timelimit / 20;

    double k = SF_INITIAL_K*2;
    double best_obj, prev_obj;
    CPXgetobjval(env, lp, &prev_obj);

    while (inst->params->timelimit > 0) {
        /* track the time for the current iteration */
        struct timespec s, e;
        s.tv_sec = e.tv_sec = -1;
        stopwatch(&s, &e);

        /* set the new timelimit and perform another resolution */
        CPXsetdblparam(env, CPX_PARAM_TILIM, iter_time_limit);

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
        if (k >= 2*SF_MAX_K) break;
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

    inst->params->timelimit = original_timelimit;
}
