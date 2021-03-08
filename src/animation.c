#include "animation.h"
#include "halfedge.h"
#include "globals.h"
#include "point.h"
#include <stdlib.h>


/* drawing helper helper */
void bov_draw_line(bov_window window, point origin, point destination,
        double scatter,
        double width, color line_color);


/* drawing helpers */
void bov_draw_beachline(bov_window window, rbtree bline, double sline_height, 
        double scatter, double init_x, double target_x, 
        double width, color bline_color) {

    /* retrive the first arc in the linked list */
    rbtree_node arc = rbtree_start_linkedlist(bline);

    /* reserve the points */
    size_t num_scatters = (size_t) ((target_x - init_x) / scatter);
    GLfloat (*bline_points)[2] = malloc(sizeof(bline_points[0])*num_scatters);

    /* fill the points */
    for (size_t i=0; i<num_scatters; i++) {
        double x = init_x + i*scatter;
        bline_points[i][0] = x;

        bline_points[i][1] = GLOBALS.MAX_VALUE;
        while (arc != NULL) {
            point focus = arc->site;
            arc = arc->next;
            bline_points[i][1] = min(bline_points[i][1], 
                    (GLfloat) evaluate_parabola(focus, sline_height, x));
        }

        arc = rbtree_start_linkedlist(bline);
    }

    /* create bov points */
    bov_points bline_pointset =  bov_points_new(bline_points, num_scatters, GL_STATIC_DRAW);

    /* set params */
    bov_points_set_color(bline_pointset, to_bov_color(bline_color));
    bov_points_set_width(bline_pointset, width);

    /* draw */
    bov_line_strip_draw(window, bline_pointset, 0, bline_pointset->vboLen);
    /* free after delation */
    bov_points_delete(bline_pointset);
}
void bov_draw_pointset(bov_window window, point* sites, size_t num_sites,
        double width, color points_color) {

    GLfloat (*sites_points)[2] = malloc(sizeof(sites_points[0])*num_sites);
    for (size_t i=0; i<num_sites; i++) {
        sites_points[i][0] = sites[i].x;
        sites_points[i][1] = sites[i].y;
    }

    /* create bov points */
    bov_points site_pointset = bov_points_new(sites_points, num_sites, GL_STATIC_DRAW);

    /* set params */
    bov_points_set_color(site_pointset, to_bov_color(points_color));
    bov_points_set_width(site_pointset, width);

    /* draw */
    bov_points_draw(window, site_pointset, 0, site_pointset->vboLen);
    /* free after delation */
    bov_points_delete(site_pointset);
}
void bov_draw_sweepline(bov_window window, double sline_height,
        double scatter, double init_x, double target_x, 
        double width, color sline_color) {

    /* reserve the points */
    size_t num_scatters = (size_t) ((target_x - init_x) / scatter);
    GLfloat (*sline_points)[2] = malloc(sizeof(sline_points[0])*num_scatters);

    /* fill the points */
    for (size_t i=0; i<num_scatters; i++) {
        double x = init_x + i*scatter;
        sline_points[i][0] = x;
        sline_points[i][1] = sline_height;
    }

    /* create bov points */
    bov_points sline_pointset = bov_points_new(sline_points, num_scatters, GL_STATIC_DRAW);

    /* set params */
    bov_points_set_color(sline_pointset, to_bov_color(sline_color));
    bov_points_set_width(sline_pointset, width);

    /* draw */
    bov_points_draw(window, sline_pointset, 0, sline_pointset->vboLen);
    /* free after delation */
    bov_points_delete(sline_pointset);
}
void bov_draw_halfedges(bov_window window, voronoi_halfedge* halfedges, size_t halfedges_size,
        double scatter,
        double width, color halfedge_color) {

    for (size_t i=0; i<halfedges_size; i++) {
        if (!halfedge_valid(halfedges[i])) continue;

        bov_draw_line(window, halfedges[i]->origin->site, halfedges[i]->destination->site,
                scatter,
                width, halfedge_color);
    }
}
void bov_draw_unvalid_halfedges(bov_window window, rbtree bline, double sline_height, 
        double scatter,
        double width, color halfedge_color) {

    /* iterate all over the linked list */
    rbtree_node act = rbtree_start_linkedlist(bline);

    /* loop n-1 times over n item in linked list */
    while (act->next != NULL) {

        /* pick first halgedge from current arc */
        voronoi_halfedge he = act->face->he;
        while (he->prev != NULL) he = he->prev;

        /* pick last halgedge from next arc */
        voronoi_halfedge next_he = act->next->face->he;
        while (next_he->next != NULL) next_he = next_he->next;

        /* pick starting vertex */
        voronoi_vertex origin = he->destination;
        voronoi_vertex origin_next = next_he->origin;


        /* if a situation of type ABA happens: draw he between parabolas */
        if (origin == NULL || origin_next == NULL) {
            point breakpoint1;
            breakpoint1.x = breakpoint_position(act->site, act->next->site, sline_height);
            breakpoint1.y = evaluate_parabola(act->site, sline_height, breakpoint1.x);

            point breakpoint2;
            breakpoint2.x = breakpoint_position(act->next->site, act->site, sline_height);
            breakpoint2.y = evaluate_parabola(act->site, sline_height, breakpoint2.x);

            bov_draw_line(window, breakpoint1, breakpoint2,
                    scatter,
                    width, halfedge_color);
        }
        /*else if (origin != origin_next) {*/
            /*double a, b, c;*/

            /*line_params(act->site, act->next->site, &a, &b, &b);*/
            /*double m_ref = -b/a;*/

            /*point breakpoint;*/
            /*breakpoint.x = breakpoint_position(act->site, act->next->site, sline_height);*/
            /*breakpoint.y = evaluate_parabola(act->site, sline_height, breakpoint.x);*/

            /*line_params(breakpoint, origin->site, &a, &b, &b);*/
            /*double m1 = -b/a;*/

            /*line_params(breakpoint, origin_next->site, &a, &b, &b);*/
            /*double m2 = -b/a;*/

            /*if (fabs(m_ref - 1/m1) < fabs(m_ref - 1/m2)){*/
                /*bov_draw_line(window, breakpoint, origin->site,*/
                        /*scatter,*/
                        /*width, halfedge_color);*/
            /*}*/
            /*else {*/
                /*bov_draw_line(window, breakpoint, origin_next->site,*/
                        /*scatter,*/
                        /*width, halfedge_color);*/
            /*}*/

        /*}*/
        /* else draw not complete halfedge */
        else if (origin == origin_next) {
            point breakpoint;
            breakpoint.x = breakpoint_position(act->site, act->next->site, sline_height);
            breakpoint.y = evaluate_parabola(act->site, sline_height, breakpoint.x);


            bov_draw_line(window, breakpoint, origin->site,
                    scatter,
                    width, halfedge_color);
        }
        act = act->next;
    }
}
void bov_draw_line(bov_window window, point origin, point destination,
        double scatter,
        double width, color line_color) {

        GLfloat (*halfedge_points)[2] = malloc(sizeof(halfedge_points[0])*2);
        halfedge_points[0][0] = origin.x;
        halfedge_points[0][1] = origin.y;
        halfedge_points[1][0] = destination.x;
        halfedge_points[1][1] = destination.y;


        bov_points halfedge_pointset = bov_points_new(halfedge_points, 2, GL_STATIC_DRAW);

        bov_points_set_color(halfedge_pointset, to_bov_color(line_color));
        bov_points_set_width(halfedge_pointset, width);

        bov_line_strip_draw(window, halfedge_pointset, 0, halfedge_pointset->vboLen);
        /* free after delation */
        bov_points_delete(halfedge_pointset);
}

/* helpers */
GLfloat* to_bov_color(color rgb_color) {
    GLfloat* bov_color = malloc(4*sizeof(GLfloat));
    bov_color[0] = rgb_color.r / 255.0;
    bov_color[1] = rgb_color.g / 255.0;
    bov_color[2] = rgb_color.b / 255.0;
    bov_color[3] = 1.0;

    return  bov_color;
}



/*bov_points_t* bov_halfedges_points(double scatter, voronoi_halfedge* halfedges, size_t halfedges_size) {*/

    /*size_t num_total_scatters = 0;*/
    /*for (size_t i=0; i<halfedges_size; i++) {*/
        /*point p1 = halfedges[i]->origin->site,*/
              /*p2 = halfedges[i]->destination->site;*/
        /*num_total_scatters += (size_t) (fabs(p2.x-p1.x) / scatter);*/
    /*}*/

    /*GLfloat (*halfedges_points)[2] = malloc(sizeof(halfedges_points[0])*num_total_scatters*10);*/
    /*size_t j = 0;*/

    /*for (size_t i=0; i<halfedges_size; i++) {*/
        /*double a, b, c;*/

        /*point p1 = halfedges[i]->origin->site,*/
              /*p2 = halfedges[i]->destination->site;*/

		/*if (fabs(p1.x-p2.x) < GLOBALS.eps) {*/
			/*a = 1.0;*/
			/*b = 0.0;*/
			/*c = -p1.x;*/
		/*}*/
		/*else {*/
			/*a = -(p1.y-p2.y)/(p1.x-p2.x);*/
			/*b = 1.0;*/
			/*c = -a*p1.x - p1.y;*/
		/*}*/

        /*size_t num_scatters = (size_t) (fabs(p2.x-p1.x) / scatter);*/

        /*double init_x = min(p1.x, p2.x),*/
               /*target_x = max(p1.x, p2.x);*/

        /*[> fill the points <]*/
        /*double x = init_x;*/
        /*size_t k = 0;*/
        /*while (x < target_x) {*/
            /*double x = init_x + k*scatter;*/

            /*halfedges_points[j][0] = x;*/
            /*halfedges_points[j][1] = -(a*x+c)/b;*/

            /*j++;*/
        /*}*/
    /*}*/

    /*[> resize the array <]*/
    /*halfedges_points = realloc(halfedges_points, sizeof(halfedges_points[0])*j);*/

    /*return bov_points_new(halfedges_points, j, GL_STATIC_DRAW);*/
/*}*/
/*bov_points_t* bov_face_points(double scatter, voronoi_face f) {*/

    /*size_t num_halfedges = 0;*/
    /*voronoi_halfedge first = f->he, act = f->he;*/

    /*printf("faces\n");*/
    /*size_t num_total_scatters = (size_t) (fabs(act->origin->site.x - first->origin->site.x) / scatter);*/
    /*do {*/
        /*printpt(act->origin->site,  NOLF);*/
        /*printf(" to ");*/
        /*printpt(act->destination->site,  LF);*/

        /*point p1 = act->origin->site,*/
              /*p2 = act->destination->site;*/
        /*num_total_scatters += (size_t) (fabs(p2.x-p1.x) / scatter);*/

        /*num_halfedges++;*/
        /*act = act->next;*/
    /*} while (act != NULL && act != first);*/
    /*printf("\n\n");*/

    /*GLfloat (*faces_points)[2] = malloc(sizeof(faces_points[0])*num_total_scatters*10);*/
    /*size_t j = 0;*/

    /*act = first;*/
    /*for (size_t i=0; i<num_halfedges; i++) {*/
        /*double a, b, c;*/

        /*point p1 = act->origin->site,*/
              /*p2 = act->destination->site;*/

		/*if (fabs(p1.x-p2.x) < GLOBALS.eps) {*/
			/*a = 1.0;*/
			/*b = 0.0;*/
			/*c = -p1.x;*/
		/*}*/
		/*else {*/
			/*a = -(p1.y-p2.y)/(p1.x-p2.x);*/
			/*b = 1.0;*/
			/*c = -a*p1.x - p1.y;*/
		/*}*/

        /*size_t num_scatters = (size_t) (fabs(p2.x-p1.x) / scatter);*/

        /*double init_x = min(p1.x, p2.x),*/
               /*target_x = max(p1.x, p2.x);*/

        /*[> fill the points <]*/
        /*double x = init_x;*/
        /*size_t k = 0;*/
        /*while (x < target_x) {*/
            /*x = init_x + k*scatter;*/

            /*faces_points[j][0] = x;*/
            /*faces_points[j][1] = -(a*x+c)/b;*/

            /*j++;*/
            /*k++;*/
        /*}*/
        /*act = act->next;*/
    /*}*/

    /*[> resize the array <]*/
    /*faces_points = realloc(faces_points, sizeof(faces_points[0])*j);*/

    /*return bov_points_new(faces_points, j, GL_STATIC_DRAW);*/
/*}*/
