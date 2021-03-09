#ifndef _EVENT_H_
#define _EVENT_H_

#include "point.h"
#include "halfedge.h"

/* forward declaration */
typedef struct rbtree_node_t* rbtree_node;

enum event_type {SITE, CIRCLE};
enum event_delete {KEEP, DELETE};

typedef struct event_t {
    enum event_type type;
    enum event_delete delete;

    point site;
    voronoi_face face;
    rbtree_node node;
} *event;

/* event utils */
event event_create(enum event_type etype, enum event_delete edelete, 
        point site, voronoi_face face, rbtree_node node);
void event_print(event ev);

#endif /* _EVENT_H_ */
