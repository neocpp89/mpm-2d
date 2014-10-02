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

void bc_init(job_t *job)
{
    return;
}

/* Returns 0 for failure, anything else for success (typically 1). */
int bc_validate(job_t *job)
{
#if 0
    if (job->boundary.num_int_props != 3) {
        fprintf(stderr, "%s:%s: Wrong number of integer boundary condition properties. Expected 3, got %d.\n", __FILE__, __func__, job->boundary.num_int_props);
        return 0;
    } else {
        printf("Using offset row %d and offset column %d.\n",
            job->boundary.int_props[0], job->boundary.int_props[1]);
    }
#endif

    if (job->boundary.num_fp64_props != 2) {
        fprintf(stderr, "%s:%s: Wrong number of floating point boundary condition properties. Expected 2, got %d.\n", __FILE__, __func__, job->boundary.num_fp64_props);
        return 0;
    }

    printf("Hole radius: %f\nOpen time: %f\n",
        job->boundary.fp64_props[0], job->boundary.fp64_props[1]);

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

    /* hole radius is given by first property, time of trapdoor by second */
    double hole_radius;
    double open_time;

    /**/
    int off_row = job->N - 1;
    int off_col = job->N - 1;
/*    int off_col = ((job->N - 1)/ 2);*/

    hole_radius = job->boundary.fp64_props[0];
    open_time = job->boundary.fp64_props[1];

/*    off_row = job->boundary.int_props[0];*/
/*    off_col = job->boundary.int_props[1];*/

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

        /* bottom of hourglass */
        job->u_dirichlet[NODAL_DOF * n + XDOF_IDX] = 0;
        job->u_dirichlet[NODAL_DOF * n + YDOF_IDX] = 0;
        job->u_dirichlet_mask[NODAL_DOF * n + XDOF_IDX] = 1;
        job->u_dirichlet_mask[NODAL_DOF * n + YDOF_IDX] = 1;

#if 0
        /* hourglass hole */
        job->u_dirichlet[NODAL_DOF * (n + (job->N / 2) * job->N) + XDOF_IDX] = 0;
        job->u_dirichlet[NODAL_DOF * (n + (job->N / 2) * job->N) + YDOF_IDX] = 0;
        job->u_dirichlet_mask[NODAL_DOF * (n + (job->N / 2) * job->N) + XDOF_IDX] = 1;
        job->u_dirichlet_mask[NODAL_DOF * (n + (job->N / 2) * job->N) + YDOF_IDX] = 1;
#endif

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

        job->u_dirichlet[NODAL_DOF * (off_col + n*job->N) + XDOF_IDX] = 0;
        job->u_dirichlet_mask[NODAL_DOF * (off_col + n*job->N) + XDOF_IDX] = 1;

/*        job->u_dirichlet[NODAL_DOF * (off_col + n*job->N) + YDOF_IDX] = 0;*/
/*        job->u_dirichlet_mask[NODAL_DOF * (off_col + n*job->N) + YDOF_IDX] = 1;*/
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

    /* handle periodicity */
    for (i = 0; i < job->num_nodes; i++) {
        for (j = 0; j < NODAL_DOF; j++) {
            if (job->node_number_override[NODAL_DOF * i + j] != NODAL_DOF * i + j) {
                /* tie nodes together. */
                m = job->node_number_override[NODAL_DOF * i + j];
                n = NODAL_DOF * i + j;

                if (n % NODAL_DOF == XDOF_IDX) {
                    if (m % NODAL_DOF == XDOF_IDX) {
                        job->nodes[(m / NODAL_DOF)].mx_t += job->nodes[(n / NODAL_DOF)].mx_t;
                        job->nodes[(n / NODAL_DOF)].mx_t = job->nodes[(m / NODAL_DOF)].mx_t;

                        job->nodes[(m / NODAL_DOF)].mx_tt += job->nodes[(n / NODAL_DOF)].mx_tt;
                        job->nodes[(n / NODAL_DOF)].mx_tt = job->nodes[(m / NODAL_DOF)].mx_tt;
                    } else if (m % NODAL_DOF == YDOF_IDX) {
                        job->nodes[(m / NODAL_DOF)].my_t += job->nodes[(n / NODAL_DOF)].mx_t;
                        job->nodes[(n / NODAL_DOF)].mx_t = job->nodes[(m / NODAL_DOF)].my_t;

                        job->nodes[(m / NODAL_DOF)].my_tt += job->nodes[(n / NODAL_DOF)].mx_tt;
                        job->nodes[(n / NODAL_DOF)].mx_tt = job->nodes[(m / NODAL_DOF)].my_tt;
                    }
                } else if (n % NODAL_DOF == YDOF_IDX) {
                    if (m % NODAL_DOF == XDOF_IDX) {
                        job->nodes[(m / NODAL_DOF)].mx_t += job->nodes[(n / NODAL_DOF)].my_t;
                        job->nodes[(n / NODAL_DOF)].my_t = job->nodes[(m / NODAL_DOF)].mx_t;

                        job->nodes[(m / NODAL_DOF)].mx_tt += job->nodes[(n / NODAL_DOF)].my_tt;
                        job->nodes[(n / NODAL_DOF)].my_tt = job->nodes[(m / NODAL_DOF)].mx_tt;
                    }
                }
            }
        }
    }

    /* fix masses as well */
    for (i = 0; i < job->num_nodes; i++) {
        j = 0;
        if (job->node_number_override[NODAL_DOF * i + j] != NODAL_DOF * i + j) {
            /* tie nodes together. */
            m = job->node_number_override[NODAL_DOF * i + j];
            n = NODAL_DOF * i + j;

            job->nodes[(m / NODAL_DOF)].m += job->nodes[(n / NODAL_DOF)].m;
            job->nodes[(n / NODAL_DOF)].m = job->nodes[(m / NODAL_DOF)].m;
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

    /* handle periodicity */
    for (i = 0; i < job->num_nodes; i++) {
        for (j = 0; j < NODAL_DOF; j++) {
            if (job->node_number_override[NODAL_DOF * i + j] != NODAL_DOF * i + j) {
                /* tie nodes together. */
                m = job->node_number_override[NODAL_DOF * i + j];
                n = NODAL_DOF * i + j;

                if (n % NODAL_DOF == XDOF_IDX) {
                    if (m % NODAL_DOF == XDOF_IDX) {
                        job->nodes[(m / NODAL_DOF)].fx += job->nodes[(n / NODAL_DOF)].fx;
                        job->nodes[(n / NODAL_DOF)].fx = job->nodes[(m / NODAL_DOF)].fx;
                    } else if (m % NODAL_DOF == YDOF_IDX) {
                        job->nodes[(m / NODAL_DOF)].fy += job->nodes[(n / NODAL_DOF)].fx;
                        job->nodes[(n / NODAL_DOF)].fx = job->nodes[(m / NODAL_DOF)].fy;
                    }
                } else if (n % NODAL_DOF == YDOF_IDX) {
                    if (m % NODAL_DOF == XDOF_IDX) {
                        job->nodes[(m / NODAL_DOF)].fx += job->nodes[(n / NODAL_DOF)].fy;
                        job->nodes[(n / NODAL_DOF)].fy = job->nodes[(m / NODAL_DOF)].fx;
                    } else if (m % NODAL_DOF == YDOF_IDX) {
                        job->nodes[(m / NODAL_DOF)].fy += job->nodes[(n / NODAL_DOF)].fy;
                        job->nodes[(n / NODAL_DOF)].fy = job->nodes[(m / NODAL_DOF)].fy;
                    }
                }
            }
        }
    }

    return;
}
