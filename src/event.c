#include "event.h"
#include "point.h"
#include <stdlib.h>
#include <stdio.h>

typedef enum event_type type;
typedef enum event_delete delete;

/* event utils */
event event_create(enum event_type etype, enum event_delete edelete, 
        point site, voronoi_face face, rbtree_node node) {

    event ev = (event) calloc(1, sizeof(struct event_t));
    ev->type = etype;
    ev->delete = edelete;
    ev->site = site;
    ev->face = face;
    ev->node = node;

    return ev;
}
void event_print(event ev) {
    char* etype_string = ev->type == SITE ? "SITE" : "CIRCLE";
    char* edelete_string = ev->delete == KEEP ? "KEEP" : "DELETE";
    char* node_string = ev->node == NULL ? "NULL" : "not NULL";

    printf("[<%s,%s>: ", etype_string, edelete_string);
    printpt(ev->site, NOLF);
    printf(" %s]\n", node_string);
}
