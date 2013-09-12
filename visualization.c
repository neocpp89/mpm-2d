/**
    \file visualization.c
    \author Sachith Dunatunga
    \date 02.06.12

    mpm_2d -- An implementation of the Material Point Method in 2D.
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "particle.h"
#include "point.h"
#include "process.h"
#include "visualization.h"

/*----------------------------------------------------------------------------*/
FILE *makeplot(int persist)
{
    FILE *pipe;

    if (persist != 0) {
        pipe = popen("gnuplot -persist", "w");
    } else {
        pipe = popen("gnuplot", "w");
    }

    fprintf(pipe, "set terminal wxt size 900,900\n");
    fprintf(pipe, "set xrange [0:1]\n");
    fprintf(pipe, "set yrange [0:1]\n");
    fprintf(pipe, "set zrange [0:1]\n");
    /* fprintf(pipe, "set cbrange [1:3]\n"); */
    fprintf(pipe, "set view map\n");
    fprintf(pipe, "set size square\n");
    fprintf(pipe, "set pointsize 0.5\n");
    fprintf(pipe, "set palette defined (0 '#000090', 1 '#000fff', ");
    fprintf(pipe, "2 '#0090ff', 3 '#0fffee', 4 '#90ff70', 5 '#ffee00', ");
    fprintf(pipe, "6 '#ff7000', 7 '#ee0000', 8 '#7f0000')\n");
/*    fprintf(pipe, "set palette defined (0 '#000000', 1 '#ffffff')\n");*/
    /* fprintf(pipe, "set palette defined (0 '#c7931e', 1 '#ffff00')\n"); */
    fprintf(pipe, "set nokey\n");
    fflush(pipe);

    return pipe;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
FILE *makeplot2(int persist, char *datafile)
{
    FILE *pipe;

    if (persist != 0) {
        pipe = popen("gnuplot -persist", "w");
    } else {
        pipe = popen("gnuplot", "w");
    }

    fprintf(pipe, "set terminal wxt size 900,900\n");
    /* fprintf(pipe, "set cbrange [1:3]\n"); */
    fprintf(pipe, "set view map\n");
    fprintf(pipe, "set size square\n");
    fprintf(pipe, "set pointsize 0.5\n");
    fprintf(pipe, "set palette defined (0 '#000090', 1 '#000fff', ");
    fprintf(pipe, "2 '#0090ff', 3 '#0fffee', 4 '#90ff70', 5 '#ffee00', ");
    fprintf(pipe, "6 '#ff7000', 7 '#ee0000', 8 '#7f0000')\n");
    /* fprintf(pipe, "set palette defined (0 '#c7931e', 1 '#ffff00')\n"); */
    fprintf(pipe, "set nokey\n");
    fflush(pipe);

    fprintf(pipe, "plot \"%s\"\n", datafile);
    fflush(pipe);

    return pipe;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
FILE *makeplot3(int persist, char *datafile)
{
    FILE *pipe;

    if (persist != 0) {
        pipe = popen("gnuplot -persist", "w");
    } else {
        pipe = popen("gnuplot", "w");
    }

    fprintf(pipe, "set terminal wxt size 900,900\n");
    /* fprintf(pipe, "set cbrange [1:3]\n"); */
    fprintf(pipe, "set view map\n");
    fprintf(pipe, "set size square\n");
    fprintf(pipe, "set pointsize 0.5\n");
    fprintf(pipe, "set palette defined (0 '#000090', 1 '#000fff', ");
    fprintf(pipe, "2 '#0090ff', 3 '#0fffee', 4 '#90ff70', 5 '#ffee00', ");
    fprintf(pipe, "6 '#ff7000', 7 '#ee0000', 8 '#7f0000')\n");
    /* fprintf(pipe, "set palette defined (0 '#c7931e', 1 '#ffff00')\n"); */
    fprintf(pipe, "set nokey\n");
    fflush(pipe);

    fprintf(pipe, "plot \"%s\" using 1:2:3 with lines palette\n", datafile);
    fflush(pipe);

    return pipe;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void plot_particles(FILE *pipe, job_t *job)
{
    int i;
    fprintf(pipe, "plot \"-\" u 1:2:3 with points palette pt 7, ");
    fprintf(pipe, "\"-\" using 1:2 with points lt 1 pt 6");
    
    for (i = 0; i < job->num_elements; i++) {
        fprintf(pipe, ", \"-\" using 1:2 with lines lt 1 linecolor rgb \"blue\"");
    }

    for (i = 0; i < job->num_particles; i++) {
        fprintf(pipe, ", \"-\" using 1:2 with lines lt 1 linecolor rgb \"black\"");
    }
    fprintf(pipe, "\n");

    for (i = 0; i < job->num_particles; i++) {
        fprintf(pipe, "%g %g %g\n", job->particles[i].x, job->particles[i].y, job->particles[i].color);
    }
    fprintf(pipe, "e\n");
    for (i = 0; i < job->num_nodes; i++) {
        fprintf(pipe, "%g %g\n", job->nodes[i].x, job->nodes[i].y);
    }
    fprintf(pipe, "e\n");
    
    for (i = 0; i < job->num_elements; i++) {
        fprintf(pipe, "%g %g\n", job->nodes[job->elements[i].nodes[0]].x, job->nodes[job->elements[i].nodes[0]].y);
        fprintf(pipe, "%g %g\n", job->nodes[job->elements[i].nodes[1]].x, job->nodes[job->elements[i].nodes[1]].y);
        fprintf(pipe, "%g %g\n", job->nodes[job->elements[i].nodes[2]].x, job->nodes[job->elements[i].nodes[2]].y);
        fprintf(pipe, "%g %g\n", job->nodes[job->elements[i].nodes[3]].x, job->nodes[job->elements[i].nodes[3]].y);
        fprintf(pipe, "e\n");
    }

    /* update the domains to draw them nicely. */
/*    update_particle_domains(job);*/
/*    update_corner_domains(job);*/

    for (i = 0; i < job->num_particles; i++) {
        /* fprintf(pipe, "%g %g\n", job->particles[i].x + job->particles[i].r1xn + job->particles[i].r2xn, job->particles[i].y + job->particles[i].r1yn + job->particles[i].r2yn);

        fprintf(pipe, "%g %g\n", job->particles[i].x + job->particles[i].r1xn - job->particles[i].r2xn, job->particles[i].y + job->particles[i].r1yn - job->particles[i].r2yn);

        fprintf(pipe, "%g %g\n", job->particles[i].x - job->particles[i].r1xn - job->particles[i].r2xn, job->particles[i].y - job->particles[i].r1yn - job->particles[i].r2yn);

        fprintf(pipe, "%g %g\n", job->particles[i].x - job->particles[i].r1xn + job->particles[i].r2xn, job->particles[i].y - job->particles[i].r1yn + job->particles[i].r2yn);

        fprintf(pipe, "%g %g\n", job->particles[i].x + job->particles[i].r1xn + job->particles[i].r2xn, job->particles[i].y + job->particles[i].r1yn + job->particles[i].r2yn);*/

/*        fprintf(pipe, "%g %g\n", job->particles[i].c1x, job->particles[i].c1y);*/
/*        fprintf(pipe, "%g %g\n", job->particles[i].c2x, job->particles[i].c2y);*/
/*        fprintf(pipe, "%g %g\n", job->particles[i].c3x, job->particles[i].c3y);*/
/*        fprintf(pipe, "%g %g\n", job->particles[i].c4x, job->particles[i].c4y);*/
/*        fprintf(pipe, "%g %g\n", job->particles[i].c1x, job->particles[i].c1y);*/
/*        fprintf(pipe, "e\n");*/
    }
    fflush(pipe);

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void plot_reload(FILE *pipe)
{

    fprintf(pipe, "replot\n");
    fflush(pipe);
    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void closeplot(FILE *pipe)
{
    fprintf(pipe, "quit\n");
    fflush(pipe);
    pclose(pipe);

    return;
}
/*----------------------------------------------------------------------------*/

