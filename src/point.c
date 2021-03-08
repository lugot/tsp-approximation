#include "point.h"
#include "globals.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/* point computation */
point intersection(point a1, point b1, point a2, point b2);


/* statically create a point */
point point_create(double x, double y) {
    point* p = (point*) malloc(sizeof(point));
    p->x = x;
    p->y = y;

    /* no points are deleated, all of them are userful, we can permit some 
     * memory leaks*/
    return *p;
}

/* print helpers */
void printpt(point p, enum printmode mode) {
	printf("(%.3lf, %.3lf)",  p.x, p.y);
	if (mode == LF) printf("\n");
}
void printpts(point *v, size_t size, enum printmode mode) {
	for(size_t i=0; i<size; i++) printpt(v[i], mode);
	printf("\n");
}

/* utils */
double dist(point a, point b){
    return hypot(a.x - b.x, a.y - b.y);
}
int same_point(point a, point b) {
    if (fabs(a.x - b.x) < GLOBALS.EPS && fabs(a.y - b.y) < GLOBALS.EPS) return 1;
    return 0;
}
double min(double a, double b) {
    if (a < b) return a;
    return b;
}
double max(double a, double b) {
    if (a > b) return a;
    return b;
}

/* point computations: circucenter, breakpoints, parabola */
point circumcenter(point pi, point pj, point pk) {
    point a1 = { 0.5*(pj.x + pk.x), 0.5*(pj.y + pk.y) },
                 a2 = { 0.5*(pi.x + pk.x), 0.5*(pi.y + pk.y) };

    point b1 = { a1.x - (pk.y - pj.y), a1.y + (pk.x - pj.x) },
                 b2 = { a2.x - (pk.y - pi.y), a2.y + (pk.x - pi.x) };

    return intersection(a1, b1, a2, b2);
}
point intersection(point a1, point b1, point a2, point b2) {

    double d = (b2.y - a2.y)*(b1.x-a1.x) - (b2.x - a2.x)*(b1.y - a1.y);
    double ua = (b2.x - a2.x)*(a1.y - a2.y) - (b2.y - a2.y)*(a1.x - a2.x);
    ua /= d;

    return (point) { a1.x + ua*(b1.x - a1.x), a1.y + ua*(b1.y - a1.y)};
}
double breakpoint_position(point site1, point site2, double bline_height) {
    double x1 = site1.x, y1 = site1.y, x2 = site2.x, y2 = site2.y;
    double d1 = 1.0/(2.0*(y1 - bline_height));
	double d2 = 1.0/(2.0*(y2 - bline_height));
	double a = d1-d2;
	double b = 2.0*(x2*d2 - x1*d1);
	double c = (y1*y1 + x1*x1 - bline_height*bline_height)*d1 - 
        (y2*y2 + x2*x2 - bline_height*bline_height)*d2;
	double delta = b*b - 4.0*a*c;

    double ans; 
    return (-b + sqrt(delta))/(2.0*a);
}
double evaluate_parabola(point focus, double directix, double x) {
    double a = focus.x,
           b = focus.y,
           c = directix;
    double num = (x-a)*(x-a) + b*b - c*c,
           denom = 2*(b-c);

    return num / denom;
}
void line_params(point site1, point site2, 
        double *a, double *b, double *c) {

    *a = -(site1.y-site2.y)/(site1.x-site2.x);
    *b = 1.0;
    *c = -*a*site1.x - site1.y;
}
