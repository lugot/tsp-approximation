#include "fortune.h"
#include "globals.h"
#include "rbtree.h"
#include "pqueue.h"
#include "point.h"
#include "event.h"
#include "halfedge.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>


/* helpers */
void remove_arc(rbtree bline, voronoi_subdivision diagram, rbtree_node arc, voronoi_vertex vertex);
rbtree_node break_arc(rbtree bline, rbtree_node arc, point site, voronoi_face face);
void check_circle_event(rbtree_node middle_arc, double sline_height, pqueue event_queue);


/* fortune implementation */
fortune_status fortune_initialize(point* sites, size_t num_points) {
    /* initialize the diagram */
    voronoi_subdivision diagram = voronoi_create_subdivision(sites, num_points);

    /* initialize the beach line */
    rbtree bline = rbtree_create();

    /* initialize the event queue */
    pqueue event_queue = pqueue_create();
    /* add events to queue */
	for(size_t i=0; i<num_points; i++){
        point site = diagram->faces[i]->site;

        event ev = event_create(SITE, KEEP, sites[i], diagram->faces[i], NULL);
        pqueue_push(event_queue, sites[i].y, ev);
	}

    fortune_status status = (fortune_status) calloc(1, sizeof(struct fortune_status_t));
    status->diagram = diagram;
    status->bline = bline;
    status->event_queue = event_queue;

    status->sline_height = pqueue_top_key(event_queue);

    return status;
}

int fortune_step(fortune_status status) {
    /* retrive the status */
    voronoi_subdivision diagram = status->diagram;
    rbtree bline = status->bline;
    pqueue event_queue = status->event_queue;

    /* if the event queue is empty we're done! */
    if (pqueue_empty(event_queue)) return 0;

    /* retrive next event */
    event ev = pqueue_pop(event_queue);
    voronoi_face face = ev->face;

    /* sweep line heigth, used as directix in parabola creation */
    double sline_height = ev->site.y;
    status->sline_height = sline_height;

    if (GLOBALS.DEBUG){ 
        printf("[DEBUG] handle event ");
        event_print(ev);
    }

    /* delete if necessary, better than delete directly from the queue */
    if (ev->delete == DELETE) {
        if (GLOBALS.DEBUG) printf("[DEBUG] delete false circle event from the queue\n");
        return 1;
    }


    /* site event handling */
    if (ev->type == SITE) {

        /* first point management */
        if (bline->root == NULL) rbtree_insert_root(bline, ev->site, face);
        /* every other point */
        else {

            /* locate the arc directly above the new site event: this will 
             * be broken in order to add the new parabola associated with the new site */
            rbtree_node above_arc = rbtree_lookup(bline, ev->site, sline_height);

            if (GLOBALS.DEBUG) {
                printf("[DEBUG] break arck at site ");
                printpt(above_arc->site, LF);
            }

            /* check for false circle events */
            if (above_arc->event != NULL) {
                if (GLOBALS.DEBUG) { 
                    printf("[DEBUG] found false circle event: ");
                    event_print(above_arc->event);
                }
                above_arc->event->delete = DELETE;
                /*above_arc->event = NULL;*/
            }

            /* broke the above arc and retrive the prev and next arcs in the beachline */
            rbtree_node middle_arc = break_arc(bline, above_arc, ev->site, face);
            rbtree_node prev_arc = middle_arc->prev;
            rbtree_node next_arc = middle_arc->next;


            /* add points to subdivision */
            prev_arc->right_halfedge = voronoi_add_halfedge(diagram, prev_arc->face);
            middle_arc->left_halfedge = voronoi_add_halfedge(diagram, middle_arc->face);
            prev_arc->right_halfedge->opposite = middle_arc->left_halfedge;
            middle_arc->left_halfedge->opposite = prev_arc->right_halfedge;

            middle_arc->right_halfedge = middle_arc->left_halfedge; 
            next_arc->left_halfedge = prev_arc->right_halfedge;
            
            /* check circle events on the right and on the left of the
             * broken arc */
            check_circle_event(prev_arc, sline_height, event_queue);
            check_circle_event(next_arc, sline_height, event_queue);
        }
    }

    /* circle event handling */
    else { /* ev.type == CIRCLE */

        /* add the vertex (center of the circle) to the subdivision */
        voronoi_vertex vertex = voronoi_add_vertex(diagram, ev->site);

        rbtree_node circle_arc = ev->node;
        rbtree_node prev_arc = circle_arc->prev;
        rbtree_node next_arc = circle_arc->next;

        /* delete events associated with prev and next arcs */
        if (prev_arc->event != NULL && prev_arc->event->type == CIRCLE) {
            prev_arc->event->delete = DELETE;
            prev_arc->event = NULL;
        }
        if (next_arc->event != NULL && next_arc->event->type == CIRCLE) {
            next_arc->event->delete = DELETE;
            next_arc->event = NULL;
        }

        /* remove arc from the beach line and create halfedges */
        remove_arc(bline, diagram, circle_arc, vertex);

        /* check circle events on the right and on the left of the
         * cirlce arc */
        check_circle_event(prev_arc, sline_height, event_queue);
        check_circle_event(next_arc, sline_height, event_queue);
    }

    /* end of event handling: debug stuffs */
    if (GLOBALS.DEBUG) {
        printf("[DEBUG] end of event queue iteration\n-tree:\n");
        rbtree_print_tree(bline);
        printf("-linked list:\n");
        rbtree_print_linkedlist(bline);
        printf("-reversed linked list:\n");
        rbtree_print_reverse_linkedlist(bline);
        printf("-event queue:\n");
        pqueue_print(event_queue);
        printf("-subdivision:\n");
        voronoi_print_faces(diagram->faces, diagram->faces_size);
        printf("----------\n\n\n");
    }

    return 1;

}

rbtree_node break_arc(rbtree bline, rbtree_node arc, point site, voronoi_face face) {

    /* retrive the next and middle arc */
    rbtree_node middle_arc = rbtree_new_empty_node(site, face);
    rbtree_node left_arc = rbtree_new_empty_node(arc->site, arc->face);
    left_arc->left_halfedge = arc->left_halfedge;
    
    /* create new empty */
    rbtree_node right_arc = rbtree_new_empty_node(arc->site, arc->face);
    right_arc->right_halfedge = arc->right_halfedge;

    /* replace the central node */
    rbtree_replace(bline, arc, middle_arc);

    /* insert before and after */
    rbtree_insert_before(bline, middle_arc, left_arc);
    rbtree_insert_after(bline, middle_arc, right_arc);

    /* free the memory*/
    free(arc);

    return middle_arc;
}

void check_circle_event(rbtree_node middle_arc, double sline_height, pqueue event_queue) {
    rbtree_node prev_arc = middle_arc->prev;
    rbtree_node next_arc = middle_arc->next;

    if (prev_arc == NULL || next_arc == NULL) return;

    /* pick circumcenter */
    point center = circumcenter(next_arc->site, middle_arc->site, prev_arc->site);
    double circle_radius = dist(center, prev_arc->site);

    /* do not create circle event if the lowest point is above the beachline */
    if (center.y - circle_radius - sline_height > GLOBALS.EPS) return;


    /* compute if the arc shrinks */
    int left_moving_right  = prev_arc->site.y < middle_arc->site.y ? 1 : 0;
    int right_moving_right = middle_arc->site.y < next_arc->site.y ? 1 : 0;
    double left_x = left_moving_right ? prev_arc->site.x : middle_arc->site.x;
    double right_x = right_moving_right ? middle_arc->site.x : next_arc->site.x;

    int arc_shrink_left = (left_x < center.x && left_moving_right) || 
        (left_x > center.x && !left_moving_right);

    int arc_shrink_rigth = (right_x < center.x && right_moving_right) || 
        (right_x > center.x && !right_moving_right);

    /* a situation of ABCBCA on C insertion: we have to check ACB and BCA can
     * produce the same center (obv) but only one of them can generate a circle
     * event */ 
    if (!arc_shrink_left ||  !arc_shrink_rigth) return;


    /* add circle event in the queue */
    event circle_event = event_create(CIRCLE, KEEP, center, NULL, middle_arc);
    middle_arc->event = circle_event;
    pqueue_push(event_queue, center.y - circle_radius, circle_event);
}
void remove_arc(rbtree bline, voronoi_subdivision diagram, rbtree_node arc, voronoi_vertex vertex) {

    /* set the destinations */
    arc->prev->right_halfedge->origin = vertex;
    arc->left_halfedge->destination = vertex;
    arc->right_halfedge->origin = vertex;
    arc->next->left_halfedge->destination = vertex;

    /* join the halfedges */
    arc->left_halfedge->next = arc->right_halfedge;
    arc->right_halfedge->prev = arc->left_halfedge;

    /* finally remove the arc from the beachline */
    rbtree_remove(bline, arc);

    /* create new edge */
    voronoi_halfedge prev_halfedge = arc->prev->right_halfedge;
    voronoi_halfedge next_halfedge = arc->next->left_halfedge;

    /* add all the two halfedges to the subdivision */
    arc->prev->right_halfedge = voronoi_add_halfedge(diagram, arc->prev->face);
    arc->next->left_halfedge = voronoi_add_halfedge(diagram, arc->next->face);
    arc->prev->right_halfedge->opposite = arc->next->left_halfedge;
    arc->next->left_halfedge->opposite = arc->prev->right_halfedge;

    /* set the origins */
    arc->prev->right_halfedge->destination = vertex;
    arc->next->left_halfedge->origin = vertex;

    /* set the previous */
    arc->prev->right_halfedge->next = prev_halfedge;
    prev_halfedge->prev = arc->prev->right_halfedge;
    next_halfedge->next = arc->next->left_halfedge;
    arc->next->left_halfedge->prev = next_halfedge;

    /* free the memory */
    free(arc);
}

