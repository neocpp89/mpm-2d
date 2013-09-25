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
#include "element.h"

#define TOL 1e-10

/*----------------------------------------------------------------------------*/
void generate_dirichlet_bcs(job_t *job)
{
    int n;
    int i, j;

    #define HOLE_RAD 0.05f

/*    int off = ((job->N - 1)/ 2);*/
    int off = job->N - 1;

    for (i = 0; i < job->num_nodes; i++) {
        for (j = 0; j < NODAL_DOF; j++) {
            job->u_dirichlet_mask[NODAL_DOF * i + j] = 0;
        }
    }

    /* Floor (and ceiling commented out). */
    for (n = 0; n < job->N; n++) {
        /* trapdoor */
/*        if (job->nodes[n].x <= HOLE_RAD && job->t > 0.5f) {*/
/*            continue;*/
/*        }*/
        job->u_dirichlet[NODAL_DOF * n + XDOF_IDX] = 0;
        job->u_dirichlet[NODAL_DOF * n + YDOF_IDX] = 0;
        job->u_dirichlet_mask[NODAL_DOF * n + XDOF_IDX] = 1;
        job->u_dirichlet_mask[NODAL_DOF * n + YDOF_IDX] = 1;

/*        job->u_dirichlet[NODAL_DOF * (job->num_nodes - n - 1) + XDOF_IDX] = 0;*/
/*        job->u_dirichlet[NODAL_DOF * (job->num_nodes - n - 1) + YDOF_IDX] = 0;*/
/*        job->u_dirichlet_mask[NODAL_DOF * (job->num_nodes - n - 1) + XDOF_IDX] = 1;*/
/*        job->u_dirichlet_mask[NODAL_DOF * (job->num_nodes - n - 1) + YDOF_IDX] = 1;*/
    }

    /* Side walls. */
    for (n = 0; n < job->N; n++) {
        job->u_dirichlet[NODAL_DOF * (n*job->N) + XDOF_IDX] = 0;
        job->u_dirichlet_mask[NODAL_DOF * (n*job->N) + XDOF_IDX] = 1;

        job->u_dirichlet[NODAL_DOF * (off + n*job->N) + XDOF_IDX] = 0;
        job->u_dirichlet_mask[NODAL_DOF * (off + n*job->N) + XDOF_IDX] = 1;
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void zero_momentum_bc(job_t *job)
{
    int n;

    #define HOLE_RAD 0.05f

    int off = ((job->N - 1)/ 2);

    for (n = 0; n < job->N; n++) {
        /* "Trapdoor" */
        if (job->nodes[n].x <= HOLE_RAD && job->t > 0.5f) {
            continue;
        }
        job->nodes[n].my_t = 0;
        job->nodes[n].mx_t = 0;
    }

    /* Side walls. */
    for (n = 0; n < job->N; n++) {
        job->nodes[n*job->N].mx_t = 0;
        job->nodes[off + n*job->N].mx_t = 0;
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void zero_force_bc(job_t *job)
{
    int n;

    int off = ((job->N - 1)/ 2);

    for (n = 0; n < job->N; n++) {
        /* "Trapdoor" */
        if (job->nodes[n].x <= HOLE_RAD && job->t > 0.5f) {
            continue;
        }
        job->nodes[n].fy = 0;
        job->nodes[n].fx = 0;
    }

    /* Side walls. */
    for (n = 0; n < job->N; n++) {
        job->nodes[n*job->N].fx = 0;
        job->nodes[off + n*job->N].fx = 0;
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
    #if 0
    int i;

    for (i = 0; i < job->num_particles; i++) {
        while (job->particles[i].x < 0) {
            job->particles[i].x += 0.5;
        }
        while (job->particles[i].x > 0.5) {
            job->particles[i].x -= 0.5;
        }
    }
    #endif

    return;
}
/*----------------------------------------------------------------------------*/


