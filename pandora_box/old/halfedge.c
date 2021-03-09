#include "halfedge.h"
#include "point.h"
#include "rbtree.h"
#include "fortune.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>


typedef voronoi_subdivision subdivision;
typedef voronoi_halfedge halfedge;
typedef voronoi_face face;
typedef voronoi_vertex vertex;


/* null object creation */
halfedge create_empty_halfedge();
vertex create_empty_vertex();

/* print helpers */
void print_face(face f);


/* object check correctness */
int looped_face(face f) {
    /* check if the face is corrected initialize */

    if (f->he == NULL) return 0;
    assert(f->he != NULL);

    halfedge first = f->he,
             act = first->next;
    
    /* loop throught the halfedges */
    while (act != NULL && act != first) act = act->next;

    if (act == first) return 1;
    return 0;
}
int halfedge_valid(halfedge he) {
    if (he->origin != NULL & he->destination != NULL) return 1;
    return 0;
}

/* null object creation */
halfedge create_empty_halfedge() {
    halfedge new_halfedge = (halfedge) calloc(1, sizeof(struct voronoi_halfedge_t));
    
    return new_halfedge;
}
vertex create_empty_vertex() {
    vertex new_vertex = (vertex) calloc(1, sizeof(struct voronoi_vertex_t));
    
    return new_vertex;
}

/* object creation utils */
subdivision voronoi_create_subdivision(point* sites, size_t num_points) {
    
    subdivision sv = (subdivision) calloc(1, sizeof(struct voronoi_subdivision_t));
    
    sv->vertices_size = 0;
    sv->vertices_capacity = 10;
    sv->vertices = (vertex*) calloc(10, sizeof(struct voronoi_vertex_t));

    sv->faces_size = num_points;
    sv->faces_capacity = num_points;
    sv->faces = (face*) calloc(num_points, sizeof(struct voronoi_face_t));

    for (size_t i=0; i<num_points; i++) {
        sv->faces[i] = (face) calloc(1, sizeof(struct voronoi_face_t));

        sv->faces[i]->site = sites[i];
        sv->faces[i]->he = NULL;
    }

    sv->halfedges_size = 0;
    sv->halfedges_capacity = 10;
    sv->halfedges = (halfedge*) calloc(10, sizeof(struct voronoi_halfedge_t));

    return sv;
}
halfedge voronoi_add_halfedge(subdivision sv, face f) {

    sv->halfedges_size++;

	if (sv->halfedges_size >= sv->halfedges_capacity) {
		sv->halfedges_capacity = 2*sv->halfedges_capacity;
		sv->halfedges = (halfedge*) realloc(sv->halfedges, 
                sv->halfedges_capacity*sizeof(struct voronoi_halfedge_t));
	}

    sv->halfedges[sv->halfedges_size-1] = create_empty_halfedge();
    sv->halfedges[sv->halfedges_size-1]->incident_face = f;
    // TODO: add iterator
    
    if (f->he == NULL) f->he = sv->halfedges[sv->halfedges_size-1];

    return sv->halfedges[sv->halfedges_size-1];
}
vertex voronoi_add_vertex(subdivision sv, point site) {

    sv->vertices_size++;

	if (sv->vertices_size >= sv->vertices_capacity) {
		sv->vertices_capacity = 2*sv->vertices_capacity;
		sv->vertices = (vertex*) realloc(sv->vertices, 
                sv->vertices_capacity*sizeof(struct voronoi_vertex_t));
	}

    sv->vertices[sv->vertices_size-1] = create_empty_vertex();
    sv->vertices[sv->vertices_size-1]->site = site;
    // TODO: add iterator
    
    return sv->vertices[sv->vertices_size-1];
}


/* print helpers */
void voronoi_print_halfedges(halfedge* halfedges, size_t halfedges_size) {

    printf("halfedges list:\n");
    for (size_t i=0; i<halfedges_size; i++) {
        printf("halfedge %zu: ", i);

        if (!halfedge_valid(halfedges[i])) {
            printf("\n");
            continue;
        }

        printpt(halfedges[i]->origin->site, NOLF);
        printf(" -> ");
        printpt(halfedges[i]->destination->site, LF);
    }

    printf("num halfedges: %zu\n", halfedges_size);
}
void voronoi_print_faces(face* faces, size_t faces_size) {

    printf("face list:\n");
    for (size_t i=0; i<faces_size; i++) {
        print_face(faces[i]);
    }

    printf("num faces: %zu\n", faces_size);
}
void print_face(face f) {

    printf("face located at ");
    printpt(f->site, NOLF);

    if (looped_face(f)) {
        printf(" (looped):\n");

        point site = f->site;

        /* check if no halfedge is present (first three points) */
        if (f->he == NULL) return;

        /* iterate over faces's halfedge like linked list */
        voronoi_halfedge first = f->he,
                         act = f->he->next; 

        if (first->origin == NULL) return;

        printpt(first->origin->site, NOLF);
        printf(" ->\n");
        while (act != first) {
            printpt(act->origin->site, NOLF);
            if (act->next != NULL) printf(" ->\n");
            act = act->next;
        }
        printpt(act->origin->site, LF);
        printf("\n");
    }
    else {
        printf(" (unlooped):\n");

        halfedge act = f->he;
        if (act == NULL) return;

        while (act->prev != NULL) act = act->prev;
        act = act->next;

        if (act == NULL) return;
        while (act->next != NULL) {
            printpt(act->origin->site, NOLF);
            if (act->next != NULL) printf(" ->\n");
            act = act->next;
        }
        printf("\n");
    }
}
void voronoi_print_vertices(vertex* vertices, size_t vertices_size) {

    printf("vertices list:\n");
    for (size_t i=0; i<vertices_size; i++) {
        printpt(vertices[i]->site, LF);
    }

    printf("num vertices: %zu\n", vertices_size);
}
void voronoi_print_results(fortune_status status) {
    subdivision diagram = status->diagram;

    voronoi_print_halfedges(diagram->halfedges, diagram->halfedges_size);
    voronoi_print_faces(diagram->faces, diagram->faces_size);
    voronoi_print_vertices(diagram->vertices, diagram->vertices_size);
}


/* link utils */
void voronoi_link_infinite_halfedges(subdivision diagram, rbtree bline, 
        double sline_height, double box_bound) {

    /* iterate all over the linked list */
    rbtree_node act = rbtree_start_linkedlist(bline);

    /* loop n-1 times over n item in linked list */
    while (act->next != NULL) {
        
        /* compute breakpoint position for direction */
        double x_breakpoint = breakpoint_position(act->site, act->next->site, sline_height);
        point breakpoint = {x_breakpoint, evaluate_parabola(act->site, sline_height, x_breakpoint)};

        /* pick first halgedge from current arc */
        halfedge he = act->face->he;
        while (he->prev != NULL) he = he->prev;

        /* pick last halgedge from next arc */
        halfedge next_he = act->next->face->he;
        while (next_he->next != NULL) next_he = next_he->next;
    

        /* pick starting vertex */
        vertex origin = he->destination;
        point origin_point = origin->site;

        /* compute vector origin->destination */
        point extern_vector = {breakpoint.x - origin->site.x,
            breakpoint.y - origin->site.y};


        /* pick tranlsation multiplication st it's outside of the box */
        double c = max((box_bound - origin_point.x)/extern_vector.x,
                (-box_bound - origin_point.x)/extern_vector.x);

        point destination_point = {origin_point.x + c*extern_vector.x,
            origin_point.y + c*extern_vector.y};


        /* link the origin and destination */
        he->origin = voronoi_add_vertex(diagram, destination_point);
        next_he->destination = he->origin;

        /* link the halfedges */
        he->opposite = next_he;
        next_he->opposite = he;

        act = act->next;
    }
}
