#ifndef _UNION_FIND_H_
#define _UNION_FIND_H_

typedef struct union_find_t {
	int *p, *rank, *size_of_set, *next;
    int num_sets;
} *union_find;

union_find uf_create(int N);
int uf_find_set(union_find uf, int i);
int uf_same_set(union_find uf, int i, int j);
int uf_set_size(union_find uf, int i);
void uf_union_set(union_find uf, int i, int j);
void uf_free(union_find uf);

#endif   /* _UNION_FIND_H_ */
