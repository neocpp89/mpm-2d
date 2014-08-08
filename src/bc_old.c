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

/*----------------------------------------------------------------------------*/
void generate_dirichlet_bcs(job_t *job)
{
    int n;
    int i, j;

    /* hole radius is given by first property, time of trapdoor by second */
    double hole_radius;
    double open_time;

    hole_radius = job->boundary.fp64_props[0];
    open_time = job->boundary.fp64_props[1];

    int off = ((job->N - 1)/ 2);
/*    int off = job->N - 1; */

    for (i = 0; i < job->num_nodes; i++) {
        for (j = 0; j < NODAL_DOF; j++) {
            job->u_dirichlet_mask[NODAL_DOF * i + j] = 0;
/*            if (j == XDOF_IDX) {*/
/*                job->u_dirichlet_mask[NODAL_DOF * i + j] = 1;*/
/*                job->u_dirichlet[NODAL_DOF * i + j] = 0; */
/*            }*/
        }
    }

    /* Floor (and ceiling commented out). */
    for (n = 0; n < job->N; n++) {
        /* trapdoor */
        if (job->nodes[n].x <= hole_radius && job->t > open_time) {
            continue;
        }

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
        job->u_dirichlet[NODAL_DOF * (0 + n*job->N) + XDOF_IDX] = 0;
        job->u_dirichlet_mask[NODAL_DOF * (0 + n*job->N) + XDOF_IDX] = 1;

/*        job->u_dirichlet[NODAL_DOF * (0 + n*job->N) + YDOF_IDX] = 0;*/
/*        job->u_dirichlet_mask[NODAL_DOF * (0+ n*job->N) + YDOF_IDX] = 1;*/

        job->u_dirichlet[NODAL_DOF * (off + n*job->N) + XDOF_IDX] = 0;
        job->u_dirichlet_mask[NODAL_DOF * (off + n*job->N) + XDOF_IDX] = 1;

/*        job->u_dirichlet[NODAL_DOF * (off + n*job->N) + YDOF_IDX] = 0;*/
/*        job->u_dirichlet_mask[NODAL_DOF * (off + n*job->N) + YDOF_IDX] = 1;*/
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void generate_node_number_override(job_t *job)
{
    int i, j;
    int off = job->N - 1;

    for (i = 0; i < job->num_nodes; i++) {
        for (j = 0; j < NODAL_DOF; j++) {
            job->node_number_override[NODAL_DOF * i + j] = (NODAL_DOF * i + j);
        }
    }

    /* x direction is periodic */
/*    for (i = 0; i < job->N; i++) {*/
/*        for (j = 0; j < NODAL_DOF; j++) {*/
/*            job->node_number_override[NODAL_DOF * (off + i*job->N) + j] =*/
/*                (NODAL_DOF * (i*job->N) + j);*/
/*        }*/
/*    }*/

    /* y direction is periodic */
/*    for (i = 0; i < job->N; i++) {*/
/*        for (j = 0; j < NODAL_DOF; j++) {*/
/*            job->node_number_override[NODAL_DOF * (job->num_nodes-job->N+i) + j] =*/
/*                (NODAL_DOF * i + j);*/
/*        }*/
/*    }*/


    return;
}
/*----------------------------------------------------------------------------*/
