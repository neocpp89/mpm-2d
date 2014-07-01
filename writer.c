/**
    \file writer.c
    \author Sachith Dunatunga
    \date 04.06.12

    mpm_2d -- An implementation of the Material Point Method in 2D.
*/
#include <stdio.h>
#include <stdlib.h>

#include "process.h"
#include "writer.h"

/*---Version 2 of output format-----------------------------------------------*/
size_t v2_write_frame(const char *directory, FILE *metafd, job_t *job,
    size_t (*particle_output_fn)(FILE *, job_t *, particle_t *),
    size_t (*element_output_fn)(FILE *, job_t *, element_t *))
{
    size_t bytes_out = 0;

    

    return bytes_out;
}

size_t v2_write_particle(FILE *fd, job_t *job, particle_t *p)
{
    size_t bytes_out = 0;

    return bytes_out;
}
/*----------------------------------------------------------------------------*/

/*---write_particle-----------------------------------------------------------*/
void inline write_particle(FILE *fd, particle_t p, double active)
{
    int i, j;
    fprintf(fd, "%f %f %f %f %f %f ",
        p.m, p.v, p.x, p.y, p.x_t, p.y_t);
    fprintf(fd, "%f %f %f %18.18lf %18.18lf %18.18lf %f %18.18lf %f",
        p.sxx, p.sxy, p.syy, p.ux, p.uy, p.state[9], p.color, p.state[10],
        (double)active);

    for (i = 0; i < 4; i++) {
        for (j = 0; j < 2; j++) {
            fprintf(fd, " %f", p.corners[i][j]);
        }
    }
    fprintf(fd, "\n");

    return;
}
/*----------------------------------------------------------------------------*/

/*---write_frame--------------------------------------------------------------*/
void write_frame(FILE *fd, int frame, double time, job_t *job)
{
    int i;

    fprintf(fd, "%d %f %d\n", frame, time, job->num_particles);
    for (i = 0; i < job->num_particles; i++) {
        write_particle(fd, job->particles[i], (double)(job->active[i]));
    }

    /* dump the entire frame to disk (or wherever) */
    fflush(fd);

    return;
}
/*----------------------------------------------------------------------------*/

/*---write_element_frame------------------------------------------------------*/
void write_element_frame(FILE *fd, int frame, double time, job_t *job)
{
    int i;
    int e;
    double *sxx_acc;
    double *sxy_acc;
    double *syy_acc;
    double *v_acc;

    double x;
    double y;

    sxx_acc = (double *)malloc(sizeof(double) * job->num_elements);
    sxy_acc = (double *)malloc(sizeof(double) * job->num_elements);
    syy_acc = (double *)malloc(sizeof(double) * job->num_elements);
    v_acc = (double *)malloc(sizeof(double) * job->num_elements);

    for (i = 0; i < job->num_elements; i++) {
        sxx_acc[i] = 0.0f;
        sxy_acc[i] = 0.0f;
        syy_acc[i] = 0.0f;
        v_acc[i] = 0.0f;
    }

    for (i = 0; i < job->num_particles; i++) {
        if (job->active[i] == 0) {
            continue;
        }
        e = job->in_element[i];
        sxx_acc[e] += job->particles[i].v * job->particles[i].sxx;
        sxy_acc[e] += job->particles[i].v * job->particles[i].sxy;
        syy_acc[e] += job->particles[i].v * job->particles[i].syy;
        v_acc[e] += job->particles[i].v;
    }

    fprintf(fd, "%d %f %d\n", frame, time, job->num_elements);

    for (i = 0; i < job->num_elements; i++) {
        node_number_to_coords(&x, &y, job->elements[i].nodes[0], job->N, job->h);
        fprintf(fd, "%f %f ", x, y);
        node_number_to_coords(&x, &y, job->elements[i].nodes[1], job->N, job->h);
        fprintf(fd, "%f %f ", x, y);
        node_number_to_coords(&x, &y, job->elements[i].nodes[2], job->N, job->h);
        fprintf(fd, "%f %f ", x, y);
        node_number_to_coords(&x, &y, job->elements[i].nodes[3], job->N, job->h);
        fprintf(fd, "%f %f ", x, y);
        fprintf(fd, "%f %f %f\n",
            (v_acc[i] > 0) ? (sxx_acc[i] / v_acc[i]) : 0.0f,
            (v_acc[i] > 0) ? (sxy_acc[i] / v_acc[i]) : 0.0f,
            (v_acc[i] > 0) ? (syy_acc[i] / v_acc[i]) : 0.0f);
    }

    free(sxx_acc);
    free(sxy_acc);
    free(syy_acc);
    free(v_acc);

    return;
}
/*----------------------------------------------------------------------------*/

/*---write_state--------------------------------------------------------------*/
void write_state(FILE *fd, job_t *job)
{
    int i, j;

    fprintf(fd, "%f %f %f\n", job->t, job->dt, job->t_stop);
    fprintf(fd, "%d %d %d\n", job->num_particles, job->num_nodes, job->num_elements);
    fprintf(fd, "%d %f\n", job->N, job->h);

    for (i = 0; i < job->num_particles; i++) {
        /* Position */
        fprintf(fd, "%f %f\n", job->particles[i].x, job->particles[i].y);

        /* Local Coordinates */
        fprintf(fd, "%f %f\n", job->particles[i].xl, job->particles[i].yl);

        /* Velocity */
        fprintf(fd, "%f %f\n", job->particles[i].x_t, job->particles[i].y_t);

        /* Mass */
        fprintf(fd, "%f\n", job->particles[i].m);

        /* Volume */
        fprintf(fd, "%f\n", job->particles[i].v);

        /* Initial volume */
        fprintf(fd, "%f\n", job->particles[i].v0);

        /* Stress */
        fprintf(fd, "%f %f %f\n",
            job->particles[i].sxx,
            job->particles[i].sxy,
            job->particles[i].syy);

        /* Stress rate */
/*        fprintf(fd, "%f %f %f\n",*/
/*            job->particles[i].sxx_t,*/
/*            job->particles[i].sxy_t,*/
/*            job->particles[i].syy_t);*/

        /* Strain rate */
        fprintf(fd, "%f %f %f %f\n",
            job->particles[i].exx_t,
            job->particles[i].exy_t,
            job->particles[i].eyy_t,
            job->particles[i].wxy_t);

        /* Jaumann stress increment */
/*        fprintf(fd, "%f %f %f\n",*/
/*            job->particles[i].dsjxx,*/
/*            job->particles[i].dsjxy,*/
/*            job->particles[i].dsjyy);*/

        /* Elastic and plastic strain increments */
/*        fprintf(fd, "%f %f %f\n",*/
/*            job->particles[i].deexx,*/
/*            job->particles[i].deexy,*/
/*            job->particles[i].deeyy);*/
/*        fprintf(fd, "%f %f %f\n",*/
/*            job->particles[i].depxx,*/
/*            job->particles[i].depxy,*/
/*            job->particles[i].depyy);*/

        /* Body forces */
        fprintf(fd, "%f %f\n", job->particles[i].bx, job->particles[i].by);

        /* Deformation gradient tensor */
        fprintf(fd, "%f %f %f %f\n",
            job->particles[i].Fxx, job->particles[i].Fxy,
            job->particles[i].Fyx, job->particles[i].Fyy);

        /* Initial particle domain vectors */
/*        fprintf(fd, "%f %f %f %f\n",*/
/*            job->particles[i].r1x0, job->particles[i].r1y0,*/
/*            job->particles[i].r2x0, job->particles[i].r2y0);*/

        /* Updated particle domain vectors */
/*        fprintf(fd, "%f %f %f %f\n",*/
/*            job->particles[i].r1xn, job->particles[i].r1yn,*/
/*            job->particles[i].r2xn, job->particles[i].r2yn);*/

        /* Corner positions */
/*        fprintf(fd, "%f %f %f %f %f %f %f %f\n",*/
/*            job->particles[i].c1x, job->particles[i].c1y,*/
/*            job->particles[i].c2x, job->particles[i].c2y,*/
/*            job->particles[i].c3x, job->particles[i].c3y,*/
/*            job->particles[i].c4x, job->particles[i].c4y);*/

        /* Displacements */
        fprintf(fd, "%f %f\n", job->particles[i].ux, job->particles[i].uy);

        /* Color used by splot visualization */
        fprintf(fd, "%f\n", job->particles[i].color);

        /* State Variables (for constitutive law) */
        fprintf(fd, "%d\n", DEPVAR);
        for (j = 0; j < DEPVAR; j++) {
            fprintf(fd, "%f\n", job->particles[i].state[j]);
        }

        /* Flag particle as active or not (used for discharge problems). */
        fprintf(fd, "%d\n", job->active[i]);
    }

    for (i = 0; i < job->num_nodes; i++) {
        fprintf(fd, "%f %f\n", job->nodes[i].x, job->nodes[i].y);
    }

    for (i = 0; i < job->num_elements; i++) {
         fprintf(fd, "%d %d %d %d\n",
            job->elements[i].nodes[0],
            job->elements[i].nodes[1],
            job->elements[i].nodes[2],
            job->elements[i].nodes[3]);
    }

    return;
}
/*----------------------------------------------------------------------------*/

