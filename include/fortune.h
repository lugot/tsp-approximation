#ifndef _FORTUNE_H_ 
#define _FORTUNE_H_

#include "halfedge.h"
#include "point.h"
#include "pqueue.h"
#include "rbtree.h"
#include <stddef.h>


typedef struct fortune_status_t {
    voronoi_subdivision diagram;
    rbtree bline;
    pqueue event_queue;
} *fortune_status;


/* fortune implementation */
fortune_status fortune_initialize(point* sites, size_t num_points);
int fortune_step(fortune_status status);

#endif /* _FORTUNE_H_ */
