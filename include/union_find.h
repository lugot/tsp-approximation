#ifndef INCLUDE_UNION_FIND_H_
#define INCLUDE_UNION_FIND_H_

typedef struct tree_node_t {
    int size;
    int capacity;
    int* sons;
} * tree_node;

typedef struct union_find_t {
    int *p, *rank, *size_of_set;
    int nsets;
    int N;

    tree_node* tns;
    int* stack;
    int si; /* stack index */
    int* visited;
} * union_find;

union_find uf_create(int N);
int uf_find_set(union_find uf, int i);
int uf_same_set(union_find uf, int i, int j);
int uf_set_size(union_find uf, int i);
void uf_union_set(union_find uf, int i, int j);
int uf_postorder(union_find uf, int* next);
void uf_free(union_find uf);

tree_node tn_create();
void tn_add_son(tree_node tn, int s);
void tn_free(tree_node tn);

#endif  // INCLUDE_UNION_FIND_H_
