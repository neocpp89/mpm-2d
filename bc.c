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

#define XY_TO_XDOF_IDX(x, y) (NODAL_DOF * (x + y*job->N) + XDOF_IDX)
#define XY_TO_YDOF_IDX(x, y) (NODAL_DOF * (x + y*job->N) + YDOF_IDX)

/* Returns 0 for failure, anything else for success (typically 1). */
int validate_bc_properties(job_t *job)
{
    return 1;
}

/*----------------------------------------------------------------------------*/
void generate_dirichlet_bcs(job_t *job)
{
    for (int i = 0; i < 4 * job->N; i++) {
        for (int j = 0; j < NODAL_DOF; j++) {
            job->u_dirichlet_mask[NODAL_DOF * i + j] = 0;
        }
    }

    /* Left side. */
    for (int n = 0; n < job->N; n++) {
        job->u_dirichlet[NODAL_DOF * (n * job->N) + XDOF_IDX] = 0;
        job->u_dirichlet[NODAL_DOF * (n * job->N) + YDOF_IDX] = 0;
        job->u_dirichlet_mask[NODAL_DOF * (n * job->N) + XDOF_IDX] = 1;
        job->u_dirichlet_mask[NODAL_DOF * (n * job->N) + YDOF_IDX] = 1;

    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void generate_node_number_override(job_t *job)
{
    int i, j;
    int off = 1; //single column

    for (i = 0; i < 4 * job->N; i++) {
        for (j = 0; j < NODAL_DOF; j++) {
            job->node_number_override[NODAL_DOF * i + j] = (NODAL_DOF * i + j);
        }
    }

    /* y direction is periodic */
    for (i = 0; i < job->N; i++) {
        for (j = 0; j < NODAL_DOF; j++) {
            job->node_number_override[NODAL_DOF * (off*job->N + i) + j] =
                (NODAL_DOF * (i) + j);
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/
