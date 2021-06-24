#include "../include/models/benders.h"

#include <concorde.h>
#include <string.h>

#include "../include/globals.h"
#include "../include/utils.h"

int CPXPUBLIC add_BENDERS_sec_callback_candidate(CPXCALLBACKCONTEXTptr context,
                                                 solution sol);
int CPXPUBLIC add_BENDERS_sec_callback_relaxation(CPXCALLBACKCONTEXTptr context,
                                                  solution sol);
int doit_fn_concorde(double cutval, int cutcount, int* cut, void* in);

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

void add_BENDERS_sec(CPXENVptr env, CPXLPptr lp, adjlist l) {
    /* add
     * sum i in V, j in V x_ij <= |V| - 1 for each V subtour of actual solution
     * to the model
     */
    int nnodes = l->N;

    char** cname = (char**)calloc(1, sizeof(char*));
    cname[0] = (char*)calloc(100, sizeof(char));

    /* iterate over all possible subtour nodes present in the solution and, for
     * a single subtour, add a SEC so no subtour will form during next
     * iterations */
    int* subtour;
    int subsize;
    int nsubtours = 0;
    while ((subtour = adjlist_get_subtour(l, &subsize)) != NULL) {
        nsubtours++;

        /* rhs of constraint is nodes in subtour -1, because of the other "one
         * edge entering and one exiting" for each node constrain, in a subset
         * of N nodes will be present N edges. Imposing N-1 fix this issue, at
         * least one edge has to exit the subset */
        const char sense = 'L';
        double rhs = (double)subsize - 1;

        /* fetch new row .. */
        int lastrow = CPXgetnumrows(env, lp);
        snprintf(cname[0], strlen(cname[0]), "benders_sec(%d)", subtour[0] + 1);
        /* .. and fix the sense and rhs */
        if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) {
            print_error("wrong CPXnewrows [benders]");
        }

        /* iterate over all possible ordered pair (symmetric model) of the
         * subtour and prevent that edge from forming*/
        for (int i = 0; i < subsize; i++) {
            for (int j = i + 1; j < subsize; j++) {
                if (CPXchgcoef(env, lp, lastrow,
                               xpos(subtour[i], subtour[j], nnodes), 1.0)) {
                    print_error("wrong CPXchgcoef [benders]");
                }
            }
        }

        free(subtour);
    }

    if (VERBOSE) printf("[VERBOSE] num subtour BENDERS %d\n", nsubtours);

    free(cname[0]);
    free(cname);
}

int CPXPUBLIC add_BENDERS_sec_callback_driver(CPXCALLBACKCONTEXTptr context,
                                              CPXLONG contextid,
                                              void* userhandle) {
    solution sol = (solution)userhandle;

    if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE) {
        return add_BENDERS_sec_callback_candidate(context, sol);
    }
    if (contextid == CPX_CALLBACKCONTEXT_RELAXATION) {
        return add_BENDERS_sec_callback_relaxation(context, sol);
    }

    print_error("cannot handle different contextids");
    return 1;
}

int CPXPUBLIC add_BENDERS_sec_callback_candidate(CPXCALLBACKCONTEXTptr context,
                                                 solution sol) {
    /* number of columns is n chooses 2 */
    int nedges = sol->nedges;
    int ncols = nedges * (nedges - 1) / 2;

    /* retreive actual solution */
    double* xstar = (double*)calloc(ncols, sizeof(double));
    double objval = CPX_INFBOUND;
    if (CPXcallbackgetcandidatepoint(context, xstar, 0, ncols - 1, &objval)) {
        print_error("CPXcallbackgetcandidatepoint error");
    }

    /* get node informations */
    int mynode = -1;
    int mythread = -1;
    double zbest;
    double incumbent = CPX_INFBOUND;
    CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODECOUNT, &mynode);
    CPXcallbackgetinfoint(context, CPXCALLBACKINFO_THREADID, &mythread);
    CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_SOL, &zbest);
    CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_SOL, &incumbent);
    if (VERBOSE && !SUPPRESS_CALLBACK) {
        printf("node information:\n");
        printf("- node index: %d\n", mynode);
        printf("- thread index: %d\n", mythread);
        printf("- zbest: %lf\n", zbest);
        printf("- incumbent: %lf\n", incumbent);
    }

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

    /* iterate over all possible subtour nodes present in the solution and, for
     * a single subtour, add a SEC so no subtour will form during next
     * iterations */
    int* subtour;
    int subsize;
    int nsubtours = 0;
    while ((subtour = adjlist_get_subtour(l, &subsize)) != NULL) {
        nsubtours++;

        /* rhs of constraint is nodes in subtour -1, because of the other "one
         * edge entering and one exiting" for each node constrain, in a subset
         * of N nodes will be present N edges. Imposing N-1 fix this issue, at
         * least one edge has to exit the subset */
        const char sense = 'L';
        double rhs = (double)subsize - 1;

        /* better to store indexes and positions of constraint's edges in
         * additional structures */
        int nnz = 0;
        const int izero = 0;
        int* index = (int*)malloc(ncols * sizeof(int));
        double* value = (double*)malloc(ncols * sizeof(double));

        /* iterate over all possible ordered pair (symmetric model) of the
         * subtour and prevent that edge from forming*/
        for (int i = 0; i < subsize; i++) {
            for (int j = i + 1; j < subsize; j++) {
                index[nnz] = xpos(subtour[i], subtour[j], nedges);
                value[nnz++] = 1.0;
            }
        }

        /* finally set the callback for rejecte the incumbent */
        if (rhs != nedges - 1) {
            if (CPXcallbackrejectcandidate(context, 1, nnz, &rhs, &sense,
                                           &izero, index, value)) {
                print_error("CPXcallbackrejectcandidate() error");
            }
        }

        free(index);
        free(value);
        free(subtour);
    }

    if (VERBOSE && !SUPPRESS_CALLBACK) {
        printf("[VERBOSE] num subtour BENDERS (callback) %d\n", nsubtours);
    }

    return 0;
}

int CPXPUBLIC add_BENDERS_sec_callback_relaxation(CPXCALLBACKCONTEXTptr context,
                                                  solution sol) {
    int nedges = sol->nedges;

    /* struct to pass info to doit_fn_callback */
    doit_fn_input data =
        (doit_fn_input)calloc(1, sizeof(struct doit_fn_input_t));
    data->nedges = nedges;
    data->context = context;

    /* number of columns is n chooses 2 */
    int ncols = nedges * (nedges - 1) / 2;

    /* get node informations */
    double* xstar = (double*)calloc(ncols, sizeof(double));
    double objval = CPX_INFBOUND;
    double epsilon = 0.1; /* set epsilon for CCcut_violated_cuts */
    if (CPXcallbackgetrelaxationpoint(context, xstar, 0, ncols - 1, &objval)) {
        print_error("CPXcallbackgetcandidatepoint error");
    }

    /* get node informations */
    int mynode = -1;
    int mythread = -1;
    double zbest;
    double incumbent = CPX_INFBOUND;
    CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODECOUNT, &mynode);
    CPXcallbackgetinfoint(context, CPXCALLBACKINFO_THREADID, &mythread);
    CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_SOL, &zbest);
    CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_SOL, &incumbent);
    if (VERBOSE && !SUPPRESS_CALLBACK) {
        printf("node information:\n");
        printf("- node index: %d\n", mynode);
        printf("- thread index: %d\n", mythread);
        printf("- zbest: %lf\n", zbest);
        printf("- incumbent: %lf\n", incumbent);
    }

    /* random discard the cut with 90% probability */
    // TODO(lugot): TUNE 90% hyperparameter
    if (mythread % 10 == 9) return 0;

    /* elist specify the vertices */
    int* elist = (int*)malloc(2 * ncols * sizeof(int));
    int loader = 0;
    for (int i = 0; i < nedges; i++) {
        for (int j = i + 1; j < nedges; j++) {
            elist[loader++] = i;
            elist[loader++] = j;
        }
    }

    /* componets infos */
    int ncomps = 0;
    int* comps = (int*)malloc(nedges * sizeof(int));
    int* compscount = (int*)malloc(nedges * sizeof(int));

    if (CCcut_connect_components(nedges, ncols, elist, xstar, &ncomps,
                                 &compscount, &comps)) {
        print_error("CCcut_connect_components error");
    }

    /* it seems it never happends */
    if (ncomps == 1) {
        if (VERBOSE) printf("[VERBOSE] add relaxation cut\n");
        if (CCcut_violated_cuts(nedges, ncols, elist, xstar, 2 - epsilon,
                                doit_fn_concorde, (void*)data))
            print_error("CCcut_violated_cuts error");
    }

    free(data);
    free(elist);
    free(xstar);
    free(comps);
    free(compscount);

    return 0;
}

int doit_fn_concorde(double cutval, int cutcount, int* cut, void* in) {
    doit_fn_input data = (doit_fn_input)in;
    double rhs = cutcount - 1.0;
    int nnz = 0;
    char sense = 'L';
    int purgeable = CPX_USECUT_FILTER;
    int local = 0;
    int izero = 0;

    double* value =
        (double*)calloc(cutcount * (cutcount - 1) / 2, sizeof(double));
    int* index = (int*)calloc(cutcount * (cutcount - 1) / 2, sizeof(int));

    // FIX: (data->sol)->nedges is not equal to the number of nodes!
    for (int i = 0; i < cutcount; i++) {
        for (int j = i + 1; j < cutcount; j++) {
            index[nnz] = xpos(cut[i], cut[j], data->nedges);
            value[nnz++] = 1.0;
        }
    }
    // TODO(magu): FIX segmentation fault: why?
    if (CPXcallbackaddusercuts(data->context, 1, nnz, &rhs, &sense, &izero,
                               index, value, &purgeable, &local)) {
        print_error("CPXcallbackaddusercuts() error");
    }
    // TODO(magu): TEST print cut

    free(index);
    free(value);
    return 0;
}
