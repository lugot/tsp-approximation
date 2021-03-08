#ifndef _POINT_H_
#define _POINT_H_

#include <stddef.h>

enum printmode {NOLF, LF};

typedef struct point_t {
    double x, y;
} point;


/* statically create a point */
point point_create(double x, double y);

/* print helpers */
void printpt(point p, enum printmode mode);
void printpts(point *v, size_t size, enum printmode mode);

/* utils */
double dist(point a, point b);
int same_point(point a, point b);
double min(double a, double b);
double max(double a, double b);

/* point computations: circucenter, breakpoints, parabola */
point circumcenter(point pi, point pj, point pk);
double breakpoint_position(point site1, point site2, double bline_height);
double evaluate_parabola(point focus, double directix, double x);
void line_params(point site1, point site2, 
        double *a, double *b, double *c);

#endif /* _POINT_H_ */
