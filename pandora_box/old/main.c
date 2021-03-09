#include "globals.h"
#include "inputs.h"
#include "fortune.h"
#include "halfedge.h"
#include "pqueue.h"
#include "rbtree.h"
#include "animation.h"
#include "halfedge.h"
#include "point.h"
#include <bits/getopt_core.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <assert.h>

int main(int argc, char *argv[]) {

    /* set globals */
    size_t num_sites = -1;
    GLOBALS.DEBUG = 0;
    int animation = 1;
    int final_only = 0;
    int print_output = 0;
    double animation_time = 3.0;
      
    int opt; 
    while((opt = getopt(argc, argv, "n:t:dafo")) != -1){
        switch(opt){
            case 'd':
                GLOBALS.DEBUG = 1;
                break;
            case 'n':
                num_sites = atoi(optarg);
                break;
            case 'a':
                animation = 0;
                break;
            case 'f':
                final_only = 1;
                break;
            case 'o':
                print_output = 1;
                break;
            case 't':
                animation_time = atof(optarg);
                break;
            case '?':  
                printf("unknown option: %c\n", opt);
                break;  
        }  
    }

    /* make sense to flags */
    if (final_only) animation = 0;
    if (animation) final_only = 0;

    /* rebalancing flags */
    if (final_only == 1) animation = 1;

    GLOBALS.EPS = 1e-9;
    GLOBALS.MIN_VALUE = -1e9;
    GLOBALS.MAX_VALUE = 1e9;

	int seed = (int) time(NULL);
    /*seed = 1608743581;*/
	srand(seed);
	/*printf("seed = %d\n", seed);*/

    /* generate random points */
    if (num_sites == -1) num_sites = 50;
	GLfloat (*sites)[2] = malloc(sizeof(sites[0])*num_sites);
	random_points(sites, num_sites);

    /* to simple points */
	point sites_points[num_sites];
	for(size_t i=0; i<num_sites; i++) sites_points[i] = point_create(sites[i][0], sites[i][1]); 
	/*printpts(sites_points, num_initial_points, LF);*/



    /* initialize the fortune computation */
    fortune_status status = fortune_initialize(sites_points, num_sites);
    rbtree bline = status->bline;
    voronoi_subdivision diagram = status->diagram;

    if (!animation) {
        clock_t begin = clock();

        /* run the algoritmh */
        while(fortune_step(status));

        clock_t end = clock();
        /* infinite halfedges NOT linked, they appears in the animation only */

        /* print and exit program */
        if (print_output) voronoi_print_results(status);
        printf("execution time: %f ms\n", (double)(end - begin)/CLOCKS_PER_SEC*1000);

        return 0;
    }


    /* define colors and stuff */
    drawing_colors animation_colors;
    animation_colors.background = (color) {  0,   0,   0}; /* nero */
    animation_colors.sites =      (color) {255, 129,   0}; /* arancio */
    animation_colors.vertices =   (color) {255, 255, 255}; /* bianco */
    animation_colors.halfedges =  (color) {194,   0,  65}; /* ciliegia */
    animation_colors.beachline =  (color) { 55, 128, 128}; /* mistero */
    animation_colors.sweepline =  (color) {  0,   0,  79}; /* brogna */

    drawing_widths animation_widths;
    animation_widths.sites =     0.01;
    animation_widths.vertices =  0.01;
    animation_widths.halfedges = 0.005;
    animation_widths.beachline = 0.007;
    animation_widths.sweepline = 0.007;


    double init_x = 1000;
    double target_x = -1000;
    GLfloat center[2] = {0.0, 0.0};

    /* compute center and other points stuff */
    for (size_t i=0; i<num_sites; i++) {
        init_x = min(init_x, sites[i][0]);
        target_x = max(target_x, sites[i][1]);

        center[0] -= sites[i][0];
        center[1] -= sites[i][1];
    }
    center[0] /= num_sites;
    center[1] /= num_sites;

    /* max box width */
    init_x -= 20;
    target_x += 20;


    /* window creation */
    bov_window window = bov_window_new(800, 800, "fortune animation");
    bov_window_enable_help(window);
    bov_window_set_color(window, to_bov_color(animation_colors.background));
    bov_window_translate(window, center);

    /* count number of iterations, yeah it's dumb :) */
    size_t number_iterations = 0;
    while (fortune_step(status)) number_iterations++;
    status = fortune_initialize(sites_points, num_sites);
    bline = status->bline;
    diagram = status->diagram;


    double scatter = 0.01;
    int infinite_halfedges_builded = 0;
    double u_refesh_period = (animation_time / number_iterations)*1000000;

    int last_iteration = 0;
    double sline_height;

    while(!bov_window_should_close(window)) {

        /* sorry, this flags are a bit verbose.. */
        if (fortune_step(status) == 0) {
            last_iteration = 1;

            if (infinite_halfedges_builded == 0) {

                /* not needed, actually, if we draw the unvalid halfedges! */
                /*voronoi_link_infinite_halfedges(diagram, status->bline, */
                        /*pqueue_top_key(status->event_queue)-1.0, 500);*/

                /* flag that the last halfedges are now builded */
                infinite_halfedges_builded = 1;

                /* permit the last animation */
                final_only = 0;
            }
        }

        /* sweep line height is next event key, that's because the y of 
         * the circle is clearly above the sweep line. For animation purpuse we
         * think it's better to keep this height so it's more nice to see */
        sline_height = last_iteration ? sline_height : pqueue_top_key(status->event_queue);

        if (final_only) continue;

        /* sweepline draw */
        bov_draw_sweepline(window, sline_height, 
                scatter, init_x, target_x, 
                animation_widths.sweepline, animation_colors.sweepline);

        /* sites draw */
        bov_draw_pointset(window, sites_points, num_sites, 
                animation_widths.sites, animation_colors.sites);

        /* halfeddge draw */
        bov_draw_halfedges(window, diagram->halfedges, diagram->halfedges_size, 
                scatter,
                animation_widths.halfedges, animation_colors.halfedges);

        /* external halfedges draw */
        bov_draw_unvalid_halfedges(window, bline, sline_height, 
                scatter,
                animation_widths.halfedges, animation_colors.halfedges);

        /* vertices draw */
        point* vertices_points = (point*) malloc(diagram->vertices_size*sizeof(point));
        for (size_t i=0; i<diagram->vertices_size; i++) {
            vertices_points[i] = diagram->vertices[i]->site;
        }
        bov_draw_pointset(window, vertices_points, diagram->vertices_size, 
                animation_widths.vertices, animation_colors.vertices);

        /* beachline draw */
        bov_draw_beachline(window, bline, sline_height, 
                scatter, init_x, target_x,
                animation_widths.beachline, animation_colors.beachline);



        usleep(u_refesh_period);

        /*bov_window_update_and_wait_events(window);*/
        bov_window_update(window);
    }

    /*bov_points_delete(diag);*/
    bov_window_delete(window);


    /* print the results */
    if (print_output) {
        voronoi_print_results(status);
    }


    printf("seed = %d\n", seed);

    return 0;
}
/* made with <3 with Dark-Powered VIm */
/* by lugot & nani */
