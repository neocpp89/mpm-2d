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

static double hole_radius = 0;
static double open_time = 0;
static size_t wall_col = 0;

void generate_dirichlet_bcs(job_t *job);
void generate_node_number_override(job_t *job);

void bc_init(job_t *job)
{
    return;
}

/* Returns 0 for failure, anything else for success (typically 1). */
int bc_validate(job_t *job)
{
    if (job->boundary.num_fp64_props != 3) {
        fprintf(stderr, "%s:%s: Wrong number of floating point boundary condition properties. Expected 3, got %d.\n", __FILE__, __func__, job->boundary.num_fp64_props);
        return 0;
    }


    hole_radius = job->boundary.fp64_props[0];
    open_time = job->boundary.fp64_props[1];
    wall_col = (int)((job->N - 1) * job->boundary.fp64_props[2]);

    printf("%s:%s: (Hole radius, Open time, Wall Column): (%g, %g, %zu)\n",
        __FILE__, __func__, hole_radius, open_time, wall_col);

    return 1;
}

void bc_time_varying(job_t *job)
{
    generate_dirichlet_bcs(job);
    generate_node_number_override(job);
    return;
}

/*----------------------------------------------------------------------------*/
void generate_dirichlet_bcs(job_t *job)
{
    int n;
    int i, j;

    for (i = 0; i < job->num_nodes; i++) {
        for (j = 0; j < NODAL_DOF; j++) {
            job->u_dirichlet_mask[NODAL_DOF * i + j] = 0;
        }
    }

    /* Floor (and ceiling commented out). */
    for (n = 0; n < job->N; n++) {

        /* trapdoor */
        if (job->nodes[n].x <= hole_radius && job->t > open_time) {
            continue;
        }

        /* bottom of silo */
        job->u_dirichlet[NODAL_DOF * n + XDOF_IDX] = 0;
        job->u_dirichlet[NODAL_DOF * n + YDOF_IDX] = 0;
        job->u_dirichlet_mask[NODAL_DOF * n + XDOF_IDX] = 1;
        job->u_dirichlet_mask[NODAL_DOF * n + YDOF_IDX] = 1;
    }

    /* Frictionless side walls. */
    for (n = 0; n < job->N; n++) {
        job->u_dirichlet[NODAL_DOF * (0 + n*job->N) + XDOF_IDX] = 0;
        job->u_dirichlet_mask[NODAL_DOF * (0 + n*job->N) + XDOF_IDX] = 1;

        job->u_dirichlet[NODAL_DOF * (wall_col + n*job->N) + XDOF_IDX] = 0;
        job->u_dirichlet_mask[NODAL_DOF * (wall_col + n*job->N) + XDOF_IDX] = 1;
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void generate_node_number_override(job_t *job)
{
    int i, j;

    for (i = 0; i < job->num_nodes; i++) {
        for (j = 0; j < NODAL_DOF; j++) {
            job->node_number_override[NODAL_DOF * i + j] = (NODAL_DOF * i + j);
        }
    }
    return;
}
/*----------------------------------------------------------------------------*/
/* Only zeros out entries for now... */
void bc_momentum(job_t *job)
{
    int i, j, m, n;

    for (i = 0; i < job->num_nodes; i++) {
        for (j = 0; j < NODAL_DOF; j++) {
            n = NODAL_DOF * i + j;
            m = job->node_number_override[n];
            if (job->u_dirichlet_mask[m] != 0) {
                /* only handle 0 displacement right now. */
                if (job->u_dirichlet[m] == 0) {
                    if (j == XDOF_IDX) {
                        job->nodes[i].mx_t = 0;
                    } else if (j == YDOF_IDX) {
                        job->nodes[i].my_t = 0;
                    }
                }
            }
        }
    }

    return;
}

/* Only zero out entries for now... */
void bc_force(job_t *job)
{
    int i, j, m, n;

    for (i = 0; i < job->num_nodes; i++) {
        for (j = 0; j < NODAL_DOF; j++) {
            n = NODAL_DOF * i + j;
            m = job->node_number_override[n];
            if (job->u_dirichlet_mask[m] != 0) {
                /* only handle 0 displacement right now. */
                if (job->u_dirichlet[m] == 0) {
                    if (j == XDOF_IDX) {
                        job->nodes[i].fx = 0;
                    } else if (j == YDOF_IDX) {
                        job->nodes[i].fy = 0;
                    }
                }
            }
        }
    }

    return;
}
