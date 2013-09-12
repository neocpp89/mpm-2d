/**
    \file bc.c
    \author Sachith Dunatunga
    \date 04.06.12

    mpm_2d -- An implementation of the Material Point Method in 2D.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "interpolate.h"
#include "particle.h"
#include "node.h"
#include "process.h"
#include "material.h"

#define TOL 1e-10

/*----------------------------------------------------------------------------*/
void zero_momentum_bc(job_t *job)
{
    int n;
    int p;
    float f;

    #define HOLE_RAD 0.05f

    int off = ((job->N - 1)/ 2);
/*    int off = job->N;*/

    for (n = 0; n < job->N; n++) {
/*        if (n >= ((job->N - 1)/ 2) - 2 && n <= ((job->N - 1)/ 2) + 2) {*/
/*            continue;*/
/*        }*/
        /* "Trapdoor" */
/*        if (job->nodes[n].x <= HOLE_RAD && job->t > 0.5f) {*/
/*            continue;*/
/*        }*/

        job->nodes[n].my_t = 0;
        job->nodes[n].mx_t = 0;
/*        job->nodes[job->num_nodes - n - 1].my_t = 0;*/
    }
    /* Side walls. */
/*    for (n = 0; n < job->N; n++) {*/
/*        job->nodes[n*job->N].mx_t = 0;*/
/*        job->nodes[off + n*job->N].mx_t = 0;*/
/*    }*/

    /* Shearing along X */
/*    for (n = 0; n < job->N; n++) {*/
/*        job->nodes[off * job->N + n].mx_t = 1e-1 * job->nodes[off * job->N + n].m;*/
/*        job->nodes[(off + 1) * job->N + n].mx_t = 1e-1 * job->nodes[(off + 1) * job->N + n].m;*/
/*    }*/

    for (n = 0; n < job->num_particles; n++) {
        if (job->t < 2.0) {
            break;
        }
        #define XVEL 1e-1

        if (job->t >= 2.0 && job->t <= 3.0) {
            f = job->t - 2.0;
        } else {
            f = 1;
        }

        if (job->particles[n].material == M_RIGID) {
            p = job->in_element[n];
#define jen(i,t) job->nodes[job->elements[p].nodes[i]].t

/*            jen(0, fx) += (f * XVEL * jen(0, m) - jen(0, mx_t)) / job->dt;*/
/*            jen(1, fx) += (f * XVEL * jen(1, m) - jen(1, mx_t)) / job->dt;*/
/*            jen(2, fx) += (f * XVEL * jen(2, m) - jen(2, mx_t)) / job->dt;*/
/*            jen(3, fx) += (f * XVEL * jen(3, m) - jen(3, mx_t)) / job->dt;*/

/*            jen(0, mx_t) = f * XVEL * jen(0, m);*/
/*            jen(1, mx_t) = f * XVEL * jen(1, m);*/
/*            jen(2, mx_t) = f * XVEL * jen(2, m);*/
/*            jen(3, mx_t) = f * XVEL * jen(3, m);*/
        }
    }

    /* X Periodic Conditions */
    for (n = 0; n < job->N; n++) {
        job->nodes[n*job->N].m += job->nodes[off + n*job->N].m;
        job->nodes[off + n*job->N].m = job->nodes[n*job->N].m;

        job->nodes[n*job->N].mx_t += job->nodes[off + n*job->N].mx_t;
        job->nodes[off + n*job->N].mx_t = job->nodes[n*job->N].mx_t;

        job->nodes[n*job->N].my_t += job->nodes[off + n*job->N].my_t;
        job->nodes[off + n*job->N].my_t = job->nodes[n*job->N].my_t;
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void zero_force_bc(job_t *job)
{
    int n;

    int off = ((job->N - 1)/ 2);
/*    int off = job->N;*/

    for (n = 0; n < job->N; n++) {
        /* if (n >= ((job->N - 1)/ 2) - 2 && n <= ((job->N - 1)/ 2) + 2) {
            continue;
        }*/
        /* "Trapdoor" */
/*        if (job->nodes[n].x <= HOLE_RAD && job->t > 0.5f) {*/
/*            continue;*/
/*        }*/
        job->nodes[n].fy = 0;
        job->nodes[n].fx = 0;
/*        job->nodes[job->num_nodes - n - 1].fy = 0;*/
    }
    /* Side walls. */
/*    for (n = 0; n < job->N; n++) {*/
/*        job->nodes[n*job->N].fx = 0;*/
/*        job->nodes[off + n*job->N].fx = 0;*/
/*    }*/

    /* Pressure along Y */
/*    for (n = 0; n < 40; n++) {*/
/*        job->particles[job->num_particles - n * 40 - 1].syy = (-1e3);*/
/*    }*/

    /* X Periodic Conditions */
    for (n = 0; n < job->N; n++) {
        job->nodes[n*job->N].fx += job->nodes[off + n*job->N].fx;
        job->nodes[off + n*job->N].fx = job->nodes[n*job->N].fx;

        job->nodes[n*job->N].fy += job->nodes[off + n*job->N].fy;
        job->nodes[off + n*job->N].fy = job->nodes[n*job->N].fy;
    }

    for (n = 0; n < job->num_particles; n++) {
        int p;
        float f;
        if (job->t < 2.0) {
            break;
        }

        if (job->t >= 2.0 && job->t <= 3.0) {
            f = job->t - 2.0;
        } else {
            f = 1;
        }

        if (job->particles[n].material == M_RIGID) {
            p = job->in_element[n];
#define jen(i,t) job->nodes[job->elements[p].nodes[i]].t

            jen(0, fx) = (f * XVEL * jen(0, m) - jen(0, mx_t)) / job->dt;
            jen(1, fx) = (f * XVEL * jen(1, m) - jen(1, mx_t)) / job->dt;
            jen(2, fx) = (f * XVEL * jen(2, m) - jen(2, mx_t)) / job->dt;
            jen(3, fx) = (f * XVEL * jen(3, m) - jen(3, mx_t)) / job->dt;

/*            printf("[%d]fx = %f\n", 0, jen(0, fx));*/
/*            printf("[%d]fx = %f\n", 1, jen(1, fx));*/
/*            printf("[%d]fx = %f\n", 2, jen(2, fx));*/
/*            printf("[%d]fx = %f\n", 3, jen(3, fx));*/

        }

        if (job->t > 3.0) {
            if (job->particles[n].material == M_RIGID) {
                p = job->in_element[n];
    #define jen(i,t) job->nodes[job->elements[p].nodes[i]].t

                jen(0, fy) = (-jen(0, my_t)) / job->dt;
                jen(1, fy) = (-jen(1, my_t)) / job->dt;
                jen(2, fy) = (-jen(2, my_t)) / job->dt;
                jen(3, fy) = (-jen(3, my_t)) / job->dt;

            }
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void move_grid(job_t *job)
{
/*    #pragma omp parallel*/
/*    {*/
        double m;
        int i;

/*        #pragma omp for*/
        for (i = 0; i < job->num_nodes; i++) {
            m = job->nodes[i].m;

            if (m > TOL) {
                job->nodes[i].x_tt = job->nodes[i].fx / m;
                job->nodes[i].y_tt = job->nodes[i].fy / m;

                job->nodes[i].mx_t += job->dt * job->nodes[i].fx;
                job->nodes[i].my_t += job->dt * job->nodes[i].fy;

                job->nodes[i].x_t = job->nodes[i].mx_t / m;
                job->nodes[i].y_t = job->nodes[i].my_t / m;
            } else {
                job->nodes[i].x_tt = 0;
                job->nodes[i].y_tt = 0;
                job->nodes[i].x_t = 0;
                job->nodes[i].y_t = 0;
            }
        }

/*    }*/
    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void correct_particle_velocities(job_t *job)
{
    #if 0

    int i, j, flag;
    int selected_elems[1600];
    int num_sel_elems;
    /* ugly hack */

    double f;

    j = 0;
    for (i = 0; i < job->num_particles; i++) {
        if (job->particles[i].material == M_RIGID) {
            selected_elems[j] = job->in_element[i];
            j++;
        }
    }
    num_sel_elems = j;

    for (i = 0; i < job->num_particles; i++) {
        flag = 0;

        for (j = 0; j < num_sel_elems; j++) {
            if (job->in_element[i] == selected_elems[j]) {
                flag = 1;
            }
        }

        if (flag && job->t >= 2.0) {

            if (job->t >= 2.0 && job->t <= 3.0) {
                f = job->t - 2.0;
            } else {
                f = 1;
            }

            job->particles[i].x_t = f * XVEL;
        }
    }

    #endif

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void correct_particle_positions(job_t *job)
{
    int i;

    for (i = 0; i < job->num_particles; i++) {
        while (job->particles[i].x < 0) {
            job->particles[i].x += 0.5;
        }
        while (job->particles[i].x > 0.5) {
            job->particles[i].x -= 0.5;
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/


