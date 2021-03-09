#ifndef _ANIMATION_H_
#define _ANIMATION_H_

#include "point.h"
#include "inputs.h"
#include "rbtree.h"


typedef struct color_t {
    double r, g, b;
} color;

typedef struct drawing_colors_t {
    struct color_t background;
    struct color_t sites;
    struct color_t vertices;
    struct color_t halfedges;
    struct color_t beachline;
    struct color_t sweepline;
} drawing_colors;

typedef struct drawing_widths_t {
    double sites;
    double vertices;
    double halfedges;
    double beachline;
    double sweepline;
} drawing_widths;


typedef bov_window_t* bov_window;
typedef bov_points_t* bov_points;


/* drawing helpers */
void bov_draw_beachline(bov_window window, rbtree bline, double sline_height, 
        double scatter, double init_x, double target_x, 
        double width, color bline_color);
void bov_draw_pointset(bov_window window, point* sites, size_t num_sites,
        double width, color points_color);
void bov_draw_sweepline(bov_window window, double sline_height,
        double scatter, double init_x, double target_x, 
        double width, color sline_color);
void bov_draw_halfedges(bov_window window, voronoi_halfedge* halfedges, size_t halfedges_size,
        double scatter,
        double width, color halfedge_color);
void bov_draw_unvalid_halfedges(bov_window window, rbtree bline, double sline_height, 
        double scatter,
        double width, color halfedge_color);

/* color helper */
GLfloat* to_bov_color(color rgb_color);


#endif /* _ANIMATION_H_ */
