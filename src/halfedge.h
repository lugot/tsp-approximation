#ifndef _HALFEDGE_H_
#define _HALFEDGE_H_

#include "point.h"
#include <stddef.h>

/* forward declarations */
typedef struct rbtree_t* rbtree;
typedef struct fortune_status_t* fortune_status;


typedef struct voronoi_vertex_t {
    point site;
} *voronoi_vertex;

typedef struct voronoi_halfedge_t {
    struct voronoi_vertex_t* origin;
    struct voronoi_vertex_t* destination;

    struct voronoi_face_t* incident_face;

    struct voronoi_halfedge_t* prev;
    struct voronoi_halfedge_t* next;
    struct voronoi_halfedge_t* opposite;
} *voronoi_halfedge;

typedef struct voronoi_face_t {
    point site;
    struct voronoi_halfedge_t* he;
} *voronoi_face;

typedef struct voronoi_subdivision_t {
    struct voronoi_vertex_t** vertices;
    struct voronoi_face_t** faces;
    struct voronoi_halfedge_t** halfedges;

    size_t vertices_size;
    size_t faces_size;
    size_t halfedges_size;

    size_t vertices_capacity;
    size_t faces_capacity;
    size_t halfedges_capacity;
} *voronoi_subdivision;


/* object check correctness */
int looped_face(voronoi_face f);
int halfedge_valid(voronoi_halfedge he);

/* object creation utils */
voronoi_subdivision voronoi_create_subdivision(point* sites, size_t num_points);
voronoi_halfedge voronoi_add_halfedge(voronoi_subdivision sv, voronoi_face f);
voronoi_vertex voronoi_add_vertex(voronoi_subdivision sv, point site);

/* print helpers */
void voronoi_print_vertices(voronoi_vertex* vertices, size_t vertices_size);
void voronoi_print_halfedges(voronoi_halfedge* halfedges, size_t halfedges_size);
void voronoi_print_faces(voronoi_face* faces, size_t faces_size);
void voronoi_print_results(fortune_status status);

/* link utils */
void voronoi_link_infinite_halfedges(voronoi_subdivision diagram, rbtree bline, 
        double sline_height, double box_bound);

#endif /* _HALFEDGE_H_ */
