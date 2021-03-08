#ifndef _RBTREE_H_
#define _RBTREE_H_

#include "point.h"
#include "event.h"
#include "halfedge.h"

enum rbtree_node_color { RED, BLACK };

typedef struct rbtree_node_t {
    point site;
    event event;

    voronoi_face face;
    voronoi_halfedge left_halfedge;
    voronoi_halfedge right_halfedge;

    struct rbtree_node_t* left;
    struct rbtree_node_t* right;
    struct rbtree_node_t* parent;

    struct rbtree_node_t* prev;
    struct rbtree_node_t* next;

    enum rbtree_node_color color;
} *rbtree_node;

typedef struct rbtree_t {
    rbtree_node root;
    size_t num_nodes;
} *rbtree;


/* node creation utils */
rbtree_node rbtree_new_empty_node(point site, voronoi_face face);

/* tree builder */
rbtree rbtree_create();

/* lookup node based on breakpoint positions */
rbtree_node rbtree_lookup(rbtree t, point site, double sline_height);

/* insertions */
void rbtree_insert_root(rbtree t, point site, voronoi_face face);
void rbtree_insert_before(rbtree t, rbtree_node n, rbtree_node m);
void rbtree_insert_after(rbtree t, rbtree_node n, rbtree_node m);
void rbtree_replace(rbtree t, rbtree_node n, rbtree_node m);
void substitute(rbtree t, rbtree_node n, rbtree_node m);

/* delations */
void rbtree_remove(rbtree t, rbtree_node n);

/* print utils */
void rbtree_print_tree(rbtree t);
void print_tree_helper(rbtree_node n, int indent);
void rbtree_print_linkedlist(rbtree t);
void rbtree_print_reverse_linkedlist(rbtree t);
rbtree_node rbtree_start_linkedlist(rbtree t);
rbtree_node rbtree_end_linkedlist(rbtree t);

#endif /* _RBTREE_H_ */
