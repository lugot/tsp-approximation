/* Original code strucutre found in Wikipedia's archives,
 * adapted to the beachline implementation.
 * http://en.literateprograms.org/Red-black_tree_(C)?oldid=19567
 *
 * It's a very robust implementation and it's very userful by the teaching POV
 * cause of userful verification method.
 *
 * A linked list is added and used for most of the arc's modification.
 * The tree is just the container which ensures O(log n) ops bounds.
 */

#include "rbtree.h"
#include "globals.h"
#include "point.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

typedef rbtree_node node;
typedef enum rbtree_node_color color;

/* kinship helpers */
static node grandparent(node n);
static node sibling(node n);
static node uncle(node n);

/* node creation utils */
node new_node(point site, voronoi_face face, color node_color, node left, node right);
node rbtree_new_empty_node(point site, voronoi_face face);

/* insertions */
void substitute(rbtree t, node n, node m);

/* delations */
node minimum_node(node n);

/* print utils */
void print_tree_helper(rbtree_node n, int indent);

/* rotations */
void rotate_left(rbtree t, node n);
void rotate_right(rbtree t, node n);

/* insertion balancing */
void insert_case1(rbtree t, node n);
void insert_case2(rbtree t, node n);
void insert_case3(rbtree t, node n);
void insert_case4(rbtree t, node n);
void insert_case5(rbtree t, node n);

/* delation balancing */
void delete_case1(rbtree t, node n);
void delete_case2(rbtree t, node n);
void delete_case3(rbtree t, node n);
void delete_case4(rbtree t, node n);
void delete_case5(rbtree t, node n);
void delete_case6(rbtree t, node n);

/* balancing properties verifiers */
void verify_linkedlist(rbtree t);
void verify_properties(rbtree t);
void verify_property_1(node n);
void verify_property_2(node root);
color node_color(node n);
void verify_property_4(node n);
void verify_property_5(node root);
void verify_property_5_helper(node n, int black_count, int* path_black_count);


/* kinship helpers */
node grandparent(node n) {
    assert (n != NULL);
    assert (n->parent != NULL); /* not the root node */
    assert (n->parent->parent != NULL); /* Not child of root */

    return n->parent->parent;
}
node sibling(node n) {
    assert (n != NULL);
    assert (n->parent != NULL); /* root node has no sibling */

    if (n == n->parent->left) return n->parent->right; 
    else                      return n->parent->left;
}
node uncle(node n) {
    assert (n != NULL);
    assert (n->parent != NULL); /* root node has no uncle */
   assert (n->parent->parent != NULL); /* children of root have no uncle */

    return sibling(n->parent);
}

/* node creation utils */
node new_node(point site, voronoi_face face, color node_color, node left, node right) {
    node result = (node) malloc(sizeof(struct rbtree_node_t));

    result->site = site;
    result->event = NULL;

    result->face = face;
    result->left_halfedge = NULL;
    result->right_halfedge = NULL;

    result->left = left;
    result->right = right;
    if (left != NULL) left->parent = result;
    if (right != NULL) right->parent = result;
    result->parent = NULL;

    result->prev = NULL;
    result->next = NULL;

    result->color = node_color;

    return result;
}
node rbtree_new_empty_node(point site, voronoi_face face) {
    return new_node(site, face, RED, NULL, NULL);
}

/* tree builder */
rbtree rbtree_create() {
    rbtree t = (rbtree) calloc(1, sizeof(struct rbtree_t));
    t->root = NULL;
    t->num_nodes = 0;
    /*verify_properties(t);*/

    return t;
}

/* lookup node based on breakpoint positions */
node rbtree_lookup(rbtree t, point site, double sline_height) {
    node act = t->root;
    int found = 0;

    /* based on breakpoint positions: the site breaks just an arc, found
     * recursively on the tree (not linked list!) part of the data structure. */
    while(!found) {
        double left_breakpoint = GLOBALS.MIN_VALUE;
        double right_breakpoint = GLOBALS.MAX_VALUE;

        if (act->prev != NULL) {
            left_breakpoint = breakpoint_position(act->prev->site, act->site, 
                    sline_height);
        }
        if (act->next != NULL) {
            right_breakpoint = breakpoint_position(act->site, act->next->site, 
                    sline_height);
        }

        int first_if = 0, second_if = 0;
        if (site.x < left_breakpoint && act->left != NULL) {
            act = act->left;
            first_if = 1;
        }
        else if (site.x > right_breakpoint && act->right != NULL) {
            act = act->right;
            second_if = 1;
        }
        else found = 1;

        /*assert(!found || site.x > left_breakpoint);*/
        /*assert(!found || site.x < right_breakpoint);*/
    }

    return act;
}

/* insertions */
void rbtree_insert_root(rbtree t, point site, voronoi_face face) {
    node inserted_node = new_node(site, face, BLACK, NULL, NULL);
    if (t->root == NULL) t->root = inserted_node;
    t->num_nodes++;

    /*verify_linkedlist(t);*/
}
void rbtree_insert_before(rbtree t, node n, node m) {
    /* insert m before n */

    if (n->left == NULL) {
        n->left = m;
        m->parent = n;
    }
    else {
        n->prev->right = m;
        m->parent = n->prev;
    }

    /* set pointers */
    m->prev = n->prev;
    if (m->prev != NULL) m->prev->next = m;
    m->next = n;
    n->prev = m;

    t->num_nodes++;

    /* perform balancing */
    /*insert_case1(t, m);*/
    /*verify_properties(t);*/
}
void rbtree_insert_after(rbtree t, node n, node m) {
    /* insert m after n */
    if (n->right == NULL) {
        n->right = m;
        m->parent = n;
    }
    else {
        n->next->left = m;
        m->parent = n->next;
    }

    /* set pointers */
    m->next = n->next;
    if (m->next != NULL) m->next->prev = m;
    m->prev = n;
    n->next = m;

    t->num_nodes++;

    /* perform balancing */
    /*insert_case1(t, m);*/
    /*verify_properties(t);*/
}
void rbtree_replace(rbtree t, node n, node m) {
    /* kindship fixing */
    substitute(t, n, m);

    /* fix childen */
    m->left = n->left;
    m->right = n->right;
    if (m->left != NULL)  m->left->parent = m;
    if (m->right != NULL) m->right->parent = m;

    /* fix linked list */
    m->prev = n->prev;
    m->next = n->next;
    if (m->prev != NULL) m->prev->next = m;
    if (m->next != NULL) m->next->prev = m;

    /* fix colors */
    m->color = n->color;

    /*verify_linkedlist(t);*/
}
void substitute(rbtree t, node n, node m) {
    if (n->parent == NULL) t->root = m;
    else if (n == n->parent->left) n->parent->left = m;
    else                           n->parent->right = m;
    if (m != NULL) m->parent = n->parent;
}

/* delations */
void rbtree_remove(rbtree t, node n) {
    node m = n;
    color m_color = m->color;
    node l;
  
    /* n has a NULL left node: subs his right with him */
    if (n->left == NULL) {
        l = n->right;
        substitute(t, n, n->right);
    }
    /* n has a NULL right node: subs his left with him */
    else if (n->right == NULL) {
        l = n->left;
        substitute(t, n, n->left);
    }
    /* n has both !NULL nodes: find next in tree strcture and substitute */
    else {
        m = minimum_node(n->right);
        m_color = m->color;
        l = m->right;

        if (m->parent == n) {
            if (l != NULL) l->parent = m;
        }
        else {
            substitute(t, m, m->right);
            if (n->right != NULL) {
                m->right = n->right;
                m->right->parent = m;
            }
        }
        substitute(t, n, m);
        m->left = n->left;
        m->left->parent = m;
        m->color = n->color;
    }

    /*if (l != NULL && node_color(l) == BLACK) delete_case1(t, l);*/

    /* fix linked list */
    if (n->prev != NULL) n->prev->next = n->next;
    if (n->next != NULL) n->next->prev = n->prev;

    t->num_nodes--;

    /*verify_properties(t);*/
}
node minimum_node(node n) {
    while (n->left != NULL) n = n->left;

    return n;
}

/* print utils */
void rbtree_print_tree(rbtree t) {
    print_tree_helper(t->root, 0);
    puts("");
}
void print_tree_helper(node n, int indent) {
    int i;
    if (n == NULL) {
        fputs("<empty tree>", stdout);
        return;
    }
    if (n->right != NULL) print_tree_helper(n->right, indent + 4);

    for(i=0; i<indent; i++) fputs(" ", stdout);
    if (n->color == BLACK) printf("%1.2lf\n", n->site.x);
    else printf("<%1.2lf>\n", n->site.x);
    if (n->left != NULL) print_tree_helper(n->left, indent + 4);
}
void rbtree_print_linkedlist(rbtree t) {
    node n = t->root;
    if (n == NULL) {
        printf("<empty list>\n");
        return;
    }

    while (n->prev != NULL) n = n->prev;

    while (n != NULL) {
        printpt(n->site, NOLF);
        if (n->next != NULL) printf("->");
        n = n->next;
    }
    printf("\n\n");
}
void rbtree_print_reverse_linkedlist(rbtree t) {
    node n = t->root;
    if (n == NULL) {
        printf("<empty list>\n");
        return;
    }

    n = rbtree_end_linkedlist(t);

    while (n != NULL) {
        printpt(n->site, NOLF);
        if (n->prev != NULL) printf("->");
        n = n->prev;
    }
    printf("\n\n");
}
node rbtree_start_linkedlist(rbtree t) {
    node n = t->root;
    while (n->prev != NULL) n = n->prev;

    return n;
}
node rbtree_end_linkedlist(rbtree t) {
    node n = t->root;
    while (n->next != NULL) n = n->next;

    return n;
}

/* rotations */
void rotate_left(rbtree t, node n) {
    node r = n->right;
    substitute(t, n, r);
    n->right = r->left;
    if (r->left != NULL) r->left->parent = n;
    r->left = n;
    n->parent = r;
}
void rotate_right(rbtree t, node n) {
    node L = n->left;
    substitute(t, n, L);
    n->left = L->right;
    if (L->right != NULL) L->right->parent = n;
    L->right = n;
    n->parent = L;
}

/* insertion balancing */
void insert_case1(rbtree t, node n) {
    if (n->parent == NULL) n->color = BLACK;
    else insert_case2(t, n);
}
void insert_case2(rbtree t, node n) {
    if (node_color(n->parent) == BLACK) return; /* Tree is still valid */
    else insert_case3(t, n);
}
void insert_case3(rbtree t, node n) {
    if (node_color(uncle(n)) == RED) {
        n->parent->color = BLACK;
        uncle(n)->color = BLACK;
        grandparent(n)->color = RED;
        insert_case1(t, grandparent(n));
    } 
    else insert_case4(t, n);
}
void insert_case4(rbtree t, node n) {
    if (n == n->parent->right && n->parent == grandparent(n)->left) {
        rotate_left(t, n->parent);
        n = n->left;
    }
    else if (n == n->parent->left && n->parent == grandparent(n)->right) {
        rotate_right(t, n->parent);
        n = n->right;
    }
    insert_case5(t, n);
}
void insert_case5(rbtree t, node n) {
    n->parent->color = BLACK;
    grandparent(n)->color = RED;
    if (n == n->parent->left && n->parent == grandparent(n)->left) 
        rotate_right(t, grandparent(n));
    else {
    assert (n == n->parent->right && n->parent == grandparent(n)->right);
    rotate_left(t, grandparent(n));
    }
}

/* delation balancing */
void delete_case1(rbtree t, node n) {
    if (n->parent == NULL) return;
    else delete_case2(t, n);
}
void delete_case2(rbtree t, node n) {
    if (node_color(sibling(n)) == RED) {
        n->parent->color = RED;
        sibling(n)->color = BLACK;
        if (n == n->parent->left) rotate_left(t, n->parent);
        else                      rotate_right(t, n->parent);
    }
    delete_case3(t, n);
}
void delete_case3(rbtree t, node n) {
    if (node_color(n->parent) == BLACK && node_color(sibling(n)) == BLACK &&
        node_color(sibling(n)->left) == BLACK && node_color(sibling(n)->right) == BLACK) {
        sibling(n)->color = RED;
        delete_case1(t, n->parent);
    }
    else delete_case4(t, n);
}
void delete_case4(rbtree t, node n) {
    if (node_color(n->parent) == RED && node_color(sibling(n)) == BLACK &&
        node_color(sibling(n)->left) == BLACK && node_color(sibling(n)->right) == BLACK) {
        sibling(n)->color = RED;
        n->parent->color = BLACK;
    }
    else delete_case5(t, n);
}
void delete_case5(rbtree t, node n) {
    if (n == n->parent->left && node_color(sibling(n)) == BLACK &&
        node_color(sibling(n)->left) == RED && node_color(sibling(n)->right) == BLACK) {
        sibling(n)->color = RED;
        sibling(n)->left->color = BLACK;
        rotate_right(t, sibling(n));
    }
    else if (n == n->parent->right && node_color(sibling(n)) == BLACK &&
             node_color(sibling(n)->right) == RED && node_color(sibling(n)->left) == BLACK) {
        sibling(n)->color = RED;
        sibling(n)->right->color = BLACK;
        rotate_left(t, sibling(n));
    }
    delete_case6(t, n);
}
void delete_case6(rbtree t, node n) {
    sibling(n)->color = node_color(n->parent);
    n->parent->color = BLACK;
    if (n == n->parent->left) {
        assert (node_color(sibling(n)->right) == RED);
        sibling(n)->right->color = BLACK;
        rotate_left(t, n->parent);
    }
    else {
        assert (node_color(sibling(n)->left) == RED);
        sibling(n)->left->color = BLACK;
        rotate_right(t, n->parent);
    }
}

/* balancing properties verifiers */
void verify_properties(rbtree t) {
    verify_linkedlist(t);
    verify_property_1(t->root);
    verify_property_2(t->root);
    /* property 3 is implicit */
    verify_property_4(t->root);
    verify_property_5(t->root);
}
void verify_property_1(node n) {
    assert(node_color(n) == RED || node_color(n) == BLACK);
    if (n == NULL) return;
    verify_property_1(n->left);
    verify_property_1(n->right);
}
void verify_property_2(node root) {
    assert(node_color(root) == BLACK);
}
color node_color(node n) {
    return n == NULL ? BLACK : n->color;
}
void verify_property_4(node n) {
    if (node_color(n) == RED) {
        assert (node_color(n->left) == BLACK);
        assert (node_color(n->right) == BLACK);
        assert (node_color(n->parent) == BLACK);
    }
    if (n == NULL) return;
    verify_property_4(n->left);
    verify_property_4(n->right);
}
void verify_property_5(node root) {
    int black_count_path = -1;
    verify_property_5_helper(root, 0, &black_count_path);
}
void verify_property_5_helper(node n, int black_count, int* path_black_count) {
    if (node_color(n) == BLACK) black_count++;
    if (n == NULL) {
        if (*path_black_count == -1) *path_black_count = black_count;
        else assert (black_count == *path_black_count);
        return;
    }

    verify_property_5_helper(n->left, black_count, path_black_count);
    verify_property_5_helper(n->right, black_count, path_black_count);
}
void verify_linkedlist(rbtree t) {
    node n = t->root;
    if (n == NULL) return;

    n = rbtree_start_linkedlist(t);
    while (n != NULL) {
        assert(n->next == NULL || n->next->prev == n);
        n = n->next;
    }

    n = rbtree_end_linkedlist(t);
    while (n != NULL) {
        assert(n->prev == NULL || n->prev->next == n);
        n = n->prev;
    }
}
