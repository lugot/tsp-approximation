#include "../include/approximations.h"

#include <assert.h>

#include "../include/globals.h"
#include "../include/union_find.h"
#include "../include/utils.h"

solution TSPminspantree(instance inst) {
    assert(inst != NULL);

    int nnodes = inst->nnodes;
    int nedges = nnodes * (nnodes - 1) / 2;

    solution sol = create_solution(inst, MST, nnodes);
    sol->distance_time = 0.0;

    /* build weighted edge list and sort it */
    wedge* wedges = (wedge*)malloc(nedges * sizeof(struct wedge_t));
    for (int i = 0; i < nnodes; i++) {
        for (int j = i + 1; j < nnodes; j++) {
            wedges[xpos(i, j, nnodes)] = (wedge){dist(i, j, inst), i, j};
        }
    }
    qsort(wedges, nedges, sizeof(struct wedge_t), wedgecmp);

    /* itearate over the edges and create the mst */
    union_find uf = uf_create(nnodes);
    for (int k = 0; k < nedges; k++) {
        int i, j;
        i = wedges[k].i;
        j = wedges[k].j;

        if (uf_same_set(uf, i, j)) continue;

        uf_union_set(uf, i, j);
    }

    /* save the solution: pre/post order visit the tree and save the nodes
     * encountered, starting from the uf root */
    sol->zstar = 0.0;
    int prev, next;
    prev = -1;
    int i = 0;
    while (uf_postorder(uf, &next)) {
        if (prev != -1) {
            sol->edges[i] = (edge){prev, next};
            sol->zstar += dist(prev, next, inst);

            i++;
        }
        prev = next;
    }
    /* do not forget to close the loop! */
    sol->edges[i] = (edge){prev, sol->edges[0].i};
    sol->zstar += dist(prev, sol->edges[0].i, inst);

    /* refine it! */
    int* succ;
    if ((succ = edges_tosucc(sol->edges, nnodes)) == NULL) {
        print_error("solver didnt produced a tour");
    }
    // sol->zstar += twoopt_refinement(inst, succ, nnodes);
    /* sol->zstar += threeopt_refinement(inst, succ, nnodes); */
    /* store the succ as usual edges array */
    for (int j = 0; j < nnodes; j++) {
        sol->edges[j] = (edge){j, succ[j]};
    }
    free(succ);

    if (EXTRA_VERBOSE) {
        /* add the spanning tree to the solution */

        /* tree has nnodes -1 edges, realloc */
        sol->nedges += sol->nedges - 1;
        sol->edges =
            (edge*)realloc(sol->edges, sol->nedges * sizeof(struct edge_t));

        int msti = 1 + sol->nedges / 2; /* index of mst edge to add */
        int* visited = (int*)calloc(nnodes, sizeof(int));
        for (int j = 0; j < nnodes; j++) {
            if (visited[j]) continue;

            /* visit the node and store all the exiting edges */
            visited[j] = 1;
            for (int k = 0; k < uf->tns[j]->size; k++) {
                if (visited[uf->tns[j]->sons[k]]) continue;

                sol->edges[msti] = (edge){j, uf->tns[j]->sons[k]};

                msti++;
            }
        }
        free(visited);

        /* yeah i could reuse visitd and flip the bit, but this is executed just
         * one.. let me do this pls */
        int* edgecolors = (int*)calloc(sol->nedges, sizeof(int));
        for (int j = 1 + sol->nedges / 2; j < sol->nedges; j++) {
            edgecolors[j] = 1;
        }

        /* plot the instance with special code version 100 */
        plot_graphviz(sol, edgecolors, 100);
        free(edgecolors);

        /* return the correct solution */
        sol->nedges = 1 + sol->nedges / 2;
        sol->edges = realloc(sol->edges, sol->nedges * sizeof(struct edge_t));
    }

    free(wedges);
    uf_free(uf);

    return sol;
}
