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
            if (j == 2) {
                job->u_dirichlet_mask[NODAL_DOF * i + j] = 1;
            }
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

