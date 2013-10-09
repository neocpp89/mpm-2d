/**
    \file process.c
    \author Sachith Dunatunga
    \date 04.06.12

    mpm_2d -- An implementation of the Material Point Method in 2D.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "interpolate.h"
#include "element.h"
#include "particle.h"
#include "node.h"
#include "process.h"
#include "material.h"
#include <suitesparse/cs.h>

/*#include <omp.h>*/
#include <signal.h>

#define TOL 1e-16

#define signum(x) ((int)((0 < x) - (x < 0)))

#define ijton(i,j,N) ((j)*N + (i))
#define N_TO_P(j,tok,i) SMEAR(j,tok,i,h)
#define DX_N_TO_P(j,tok,c,i) (c) * SMEAR(j,tok,i,b1)
#define DY_N_TO_P(j,tok,c,i) (c) * SMEAR(j,tok,i,b2)

#define IS_VALID_ELEMENT_COORD4(r,c,N) \
    (((r) < (N)) && ((r) >= 0) && ((c) < (N)) && ((c) >= 0))

#define FILL_ELEMENT_NEIGHBOR(en,r,c,N) \
    do { \
        if (IS_VALID_ELEMENT_COORD((r),(c),(N))) { \
            en = (r)*(N) + (c); \
        } else { \
            en = -1; \
        } \
    } while(0)

#define SMEAR SMEAR4
#define ACCUMULATE ACCUMULATE4
#define ACCUMULATE_WITH_MUL ACCUMULATE_WITH_MUL4
#define WHICH_ELEMENT WHICH_ELEMENT4
#define IS_VALID_ELEMENT_COORD IS_VALID_ELEMENT_COORD4

#define __E(j,p) j->in_element[p]
#define __N(j,p,n) j->elements[__E(j,p)].nodes[n]
#define __NE(j,e,n) j->elements[e].nodes[n]

/*#define WHICH_ELEMENT4(xp,yp,N,h) \*/
/*    (floor(xp/h) + floor(yp/h)*(N-1))*/

/* XXX: ugly, fix soon */
#define WHICH_ELEMENT4(xp,yp,N,h) \
    (((xp)<=1.0 && (xp)>=0.0 && (yp)<=1.0 && (yp)>=0.0)?((floor((xp)/(h)) + floor((yp)/(h))*((N)-1))):(-1))

#define ACCUMULATE4(acc_tok,j,tok,i,n,s) \
    j->nodes[__N(j,i,0)].acc_tok += j->s ## 1[i] * j->particles[i].tok; \
    j->nodes[__N(j,i,1)].acc_tok += j->s ## 2[i] * j->particles[i].tok; \
    j->nodes[__N(j,i,2)].acc_tok += j->s ## 3[i] * j->particles[i].tok; \
    j->nodes[__N(j,i,3)].acc_tok += j->s ## 4[i] * j->particles[i].tok;

#define ACCUMULATE_WITH_MUL4(acc_tok,j,tok,i,n,s,c) \
    j->nodes[__N(j,i,0)].acc_tok += j->s ## 1[i] * j->particles[i].tok * (c); \
    j->nodes[__N(j,i,1)].acc_tok += j->s ## 2[i] * j->particles[i].tok * (c); \
    j->nodes[__N(j,i,2)].acc_tok += j->s ## 3[i] * j->particles[i].tok * (c); \
    j->nodes[__N(j,i,3)].acc_tok += j->s ## 4[i] * j->particles[i].tok * (c);

#define SMEAR4(j,tok,i,s) ( \
    j->s ## 1[i] * j->nodes[__N(j,i,0)].tok + \
    j->s ## 2[i] * j->nodes[__N(j,i,1)].tok + \
    j->s ## 3[i] * j->nodes[__N(j,i,2)].tok + \
    j->s ## 4[i] * j->nodes[__N(j,i,3)].tok \
)

#define CHECK_ACTIVE(j,i) if (j->particles[i].active == 0) { continue; }

/* Lapack function, double precision direct solver for general matrix */
extern void dgesv_(const int *N, const int *nrhs, double *A, const int *lda, int
    *ipiv, double *b, const int *ldb, int *info);

/* Lapack function, double precision direct solver for symmetric matrix */
extern void dsysv_(char* uplo, int* n, int* nrhs, double* a, int* lda,
    int* ipiv, double* b, int* ldb, double* work, int* lwork, int* info);

/* Lapack function, double precision 2-norm of vector */
extern double dnrm2_(int *n, double *x, int *incx);

/*----------------------------------------------------------------------------*/
job_t *mpm_init(int N, double h, particle_t *particles, int num_particles, double t)
{
    int i;
    int n;

    int r,c;

    job_t *job;

    job = (job_t *)malloc(sizeof(job_t));

    job->t = 0;
    job->t_stop = t;

    job->N = N;
    job->h = h;
    job->num_nodes = N*N;
    job->num_particles = num_particles;

    job->num_elements = (N - 1) * (N - 1);

    /* Copy particles from given ICs. */
    fprintf(stderr, "Each particle is %d bytes.\n", sizeof(particle_t));
    job->particles = (particle_t *)malloc(num_particles * sizeof(particle_t));
    fprintf(stderr, "%d bytes (%.2g MB) allocated for %d particles.\n",
        num_particles * sizeof(particle_t),
        num_particles * sizeof(particle_t) / 1048576.0,
        num_particles);
    memcpy(job->particles, particles, num_particles * sizeof(particle_t));

    /* Set stress, strain to zero. */
    for (i = 0; i < job->num_particles; i++) {
        job->particles[i].sxx = 0;
        job->particles[i].sxy = 0;
        job->particles[i].syy = 0;

        job->particles[i].exx_t = 0;
        job->particles[i].exy_t = 0;
        job->particles[i].eyy_t = 0;
        job->particles[i].wxy_t = 0;

        job->particles[i].Fxx = 1;
        job->particles[i].Fxy = 0;
        job->particles[i].Fyx = 0;
        job->particles[i].Fyy = 1;

        job->particles[i].ux = 0;
        job->particles[i].uy = 0;
    }
    fprintf(stderr, "Done setting initial particle data.\n");

    /* Get node coordinates. */
    fprintf(stderr, "Each node is %d bytes.\n", sizeof(node_t));
    job->nodes = (node_t *)malloc(job->num_nodes * sizeof(node_t));
    fprintf(stderr, "%d bytes (%.2g MB) allocated for %d nodes.\n",
        job->num_nodes * sizeof(node_t),
        job->num_nodes * sizeof(node_t) / 1048576.0,
        job->num_nodes);
    for (i = 0; i < job->num_nodes; i++) {
        node_number_to_coords(&(job->nodes[i].x), &(job->nodes[i].y), i, N, h);
        /* printf("preprocessing: xn=%g yn=%g\n", job->nodes[i].x, job->nodes[i].y); */
    }

    /* Get node numbering for elements. */
    fprintf(stderr, "Each element is %d bytes.\n", sizeof(element_t));
    job->elements = (element_t *)malloc(job->num_elements * sizeof(element_t));
    fprintf(stderr, "%d bytes (%.2g MB) allocated for %d elements.\n",
        job->num_elements * sizeof(element_t),
        job->num_elements * sizeof(element_t) / 1048576.0,
        job->num_elements);
    for (i = 0; i < job->num_elements; i++) {
        n = ijton(i % (job->N - 1), i / (job->N - 1), job->N);

        job->elements[i].nodes[0] = n + ijton(0,0,job->N);
        job->elements[i].nodes[1] = n + ijton(1,0,job->N);
        job->elements[i].nodes[2] = n + ijton(1,1,job->N);
        job->elements[i].nodes[3] = n + ijton(0,1,job->N);
    }

    /* Find element neighbors -- cartesian grid only! */
    for (i = 0; i < job->num_elements; i++) {
        r = i / (job->N - 1);
        c = i % (job->N - 1);

        /*
            Fill neighbor array starting with element to the right and
           moving in a postive direction (ccw).
        */
        FILL_ELEMENT_NEIGHBOR(job->elements[i].neighbors[0], r, c+1, (job->N - 1));
        FILL_ELEMENT_NEIGHBOR(job->elements[i].neighbors[1], r+1, c+1, (job->N - 1));
        FILL_ELEMENT_NEIGHBOR(job->elements[i].neighbors[2], r+1, c, (job->N - 1));
        FILL_ELEMENT_NEIGHBOR(job->elements[i].neighbors[3], r+1, c-1, (job->N - 1));
        FILL_ELEMENT_NEIGHBOR(job->elements[i].neighbors[4], r, c-1, (job->N - 1));
        FILL_ELEMENT_NEIGHBOR(job->elements[i].neighbors[5], r-1, c-1, (job->N - 1));
        FILL_ELEMENT_NEIGHBOR(job->elements[i].neighbors[6], r-1, c, (job->N - 1));
        FILL_ELEMENT_NEIGHBOR(job->elements[i].neighbors[7], r-1, c+1, (job->N - 1));

    }

    /* Allocate space for tracking element->particle map. */
    job->in_element =  (int *)malloc(job->num_particles * sizeof(int));

    /* Allocate space for interpolation functions. */
    job->h1 = (double *)malloc(job->num_particles * sizeof(double));
    job->h2 = (double *)malloc(job->num_particles * sizeof(double));
    job->h3 = (double *)malloc(job->num_particles * sizeof(double));
    job->h4 = (double *)malloc(job->num_particles * sizeof(double));
/*    job->h5 = (double *)malloc(job->num_particles * sizeof(double));*/
/*    job->h6 = (double *)malloc(job->num_particles * sizeof(double));*/
/*    job->h7 = (double *)malloc(job->num_particles * sizeof(double));*/
/*    job->h8 = (double *)malloc(job->num_particles * sizeof(double));*/
/*    job->h9 = (double *)malloc(job->num_particles * sizeof(double));*/

    job->b11 = (double *)malloc(job->num_particles * sizeof(double));
    job->b12 = (double *)malloc(job->num_particles * sizeof(double));
    job->b13 = (double *)malloc(job->num_particles * sizeof(double));
    job->b14 = (double *)malloc(job->num_particles * sizeof(double));
/*    job->b15 = (double *)malloc(job->num_particles * sizeof(double));*/
/*    job->b16 = (double *)malloc(job->num_particles * sizeof(double));*/
/*    job->b17 = (double *)malloc(job->num_particles * sizeof(double));*/
/*    job->b18 = (double *)malloc(job->num_particles * sizeof(double));*/
/*    job->b19 = (double *)malloc(job->num_particles * sizeof(double));*/

    job->b21 = (double *)malloc(job->num_particles * sizeof(double));
    job->b22 = (double *)malloc(job->num_particles * sizeof(double));
    job->b23 = (double *)malloc(job->num_particles * sizeof(double));
    job->b24 = (double *)malloc(job->num_particles * sizeof(double));
/*    job->b25 = (double *)malloc(job->num_particles * sizeof(double));*/
/*    job->b26 = (double *)malloc(job->num_particles * sizeof(double));*/
/*    job->b27 = (double *)malloc(job->num_particles * sizeof(double));*/
/*    job->b28 = (double *)malloc(job->num_particles * sizeof(double));*/
/*    job->b29 = (double *)malloc(job->num_particles * sizeof(double));*/

/*    job->phi = (double *)malloc(job->num_particles * job->num_nodes * sizeof(double));*/

    /* max size of u_grid is NODAL_DOF * number of nodes. */
    job->vec_len = NODAL_DOF * job->num_nodes;
    job->u_grid = (double *)malloc(NODAL_DOF * sizeof(double) * job->num_nodes);
    job->du_grid = (double *)malloc(NODAL_DOF * sizeof(double) * job->num_nodes);
    job->v_grid = (double *)malloc(NODAL_DOF * sizeof(double) * job->num_nodes);
    job->a_grid = (double *)malloc(NODAL_DOF * sizeof(double) * job->num_nodes);
    job->m_grid = (double *)malloc(NODAL_DOF * sizeof(double) * job->num_nodes);
    job->f_ext_grid = (double *)malloc(NODAL_DOF * sizeof(double) * job->num_nodes);
    job->f_int_grid = (double *)malloc(NODAL_DOF * sizeof(double) * job->num_nodes);
    job->q_grid = (double *)malloc(NODAL_DOF * sizeof(double) * job->num_nodes);
    job->node_u_map = (int *)malloc(NODAL_DOF * sizeof(int) * job->num_nodes);
    job->inv_node_u_map = (int *)malloc(NODAL_DOF * sizeof(int) * job->num_nodes);

    /* for dirichlet BCs */
    job->u_dirichlet = (double *)malloc(NODAL_DOF * sizeof(double) * job->num_nodes);
    job->u_dirichlet_mask = (int *)malloc(NODAL_DOF * sizeof(int) * job->num_nodes);

    /* for periodic BCs */
    job->node_number_override = (int *)malloc(job->vec_len * sizeof(int));

/*    job->kku_grid = (double *)malloc(*/
/*                        (NODAL_DOF * job->num_nodes) **/
/*                        (NODAL_DOF * sizeof(double) * job->num_nodes));*/

    for (i = 0; i < job->vec_len; i++) {
        job->u_grid[i] = 0;
        job->du_grid[i] = 0;
        job->v_grid[i] = 0;
        job->a_grid[i] = 0;
        job->m_grid[i] = 0;
        job->f_ext_grid[i] = 0;
        job->q_grid[i] = 0;
        job->f_int_grid[i] = 0;
    }
/*    for (i = 0; i < (job->vec_len * job->vec_len); i++) {*/
/*        job->kku_grid[i] = 0;*/
/*    }*/

    for (i = 0; i < job->num_particles; i++) {
        job->in_element[i] = WHICH_ELEMENT(
            job->particles[i].x, job->particles[i].y, job->N, job->h);
        if (job->in_element[i] < 0 || job->in_element[i] > job->num_elements) {
            job->particles[i].active = 0;
        } else {
            job->particles[i].active = 1;
        }
        if ((job->in_element[i] / (job->N - 1)) % 2 == 0) {
            job->particles[i].color = 0;
        } else {
            job->particles[i].color = 1;
        }
    }

    /* set particle domains */
    for (i = 0; i < job->num_particles; i++) {
        /* seems backwards, but only because loader contains
        current particle volume only. */
        job->particles[i].v0 = job->particles[i].v;
        job->particles[i].r1_initial[0] = 0.5 * sqrt(job->particles[i].v0);
        job->particles[i].r1_initial[1] = 0;
        job->particles[i].r2_initial[0] = 0;
        job->particles[i].r2_initial[1] = 0.5 * sqrt(job->particles[i].v0);
    }

    /* intialize state variables */
    material_init(job);

    /* Set default timestep. */
    job->dt = 10 * job->h * sqrt(job->particles[0].m/(job->particles[0].v * EMOD));

    /* By default, don't output energy data. */
/*    job->ke_data = NULL;*/

    /* make sure corners know where the are on first step */
    update_particle_vectors(job);
    update_corner_positions(job);

    job->use_cpdi = 0;

    return job;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void mpm_step(job_t *job)
{
    int i, j;
    stask_t task;

    /* Clear grid quantites. */
    for (i = 0; i < job->num_nodes; i++) {
        job->nodes[i].m = 0;
        job->nodes[i].mx_t = 0;
        job->nodes[i].my_t = 0;
        job->nodes[i].mx_tt = 0;
        job->nodes[i].my_tt = 0;
        job->nodes[i].x_t = 0;
        job->nodes[i].y_t = 0;
        job->nodes[i].x_tt = 0;
        job->nodes[i].y_tt = 0;
        job->nodes[i].fx = 0;
        job->nodes[i].fy = 0;

        for (j = 0; j < NODAL_DOF; j++) {
            job->u_grid[NODAL_DOF * i + j] = 0;
            job->du_grid[NODAL_DOF * i + j] = 0;
            job->f_int_grid[NODAL_DOF * i + j] = 0;
            job->f_ext_grid[NODAL_DOF * i + j] = 0;
        }
    }

    for (i = 0; i < job->num_elements; i++) {
        job->elements[i].filled = 0;
        job->elements[i].n = 0;
        job->elements[i].m = 0;
    }

    /* Figure out which element each material point is in. */
    create_particle_to_element_map(job);

    /* Calculate shape and gradient of shape functions. */
    calculate_shapefunctions(job);

    /* Map particle state to grid quantites. */
    map_to_grid(job);

    /* Calculate node velocity. */
    calculate_node_velocity(job);

    /* matrix solve */
    implicit_solve(job);

    /* Update particle position and velocity. */
    move_particles(job);

    /* update deformation gradient, vectors and corners. */
    update_deformation_gradient(job);
    update_particle_vectors(job);
    update_corner_positions(job);

    /* update volume */
    update_particle_densities(job);

    /* Increment time. */
    job->t += job->dt;

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void create_particle_to_element_map(job_t *job)
{
    int i, p;

    /* Clear all elements */
    for (i = 0; i < job->num_elements; i++) {
        job->elements[i].filled = 0;
        job->elements[i].n = 0;
        job->elements[i].m = 0;
    }

    for (i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);
        p = WHICH_ELEMENT(
            job->particles[i].x, job->particles[i].y, job->N, job->h);

        if (p != job->in_element[i]) {
            fprintf(stderr, 
                "Particle %d left element %d, now in element %d.\n", i,
                job->in_element[i], p);
        }

        /* XXX ugly hack for domain size... */
        if (job->particles[i].x < 0.0 || job->particles[i].x > 1.0
            || job->particles[i].y < 0.0 || job->particles[i].y > 1.0) {
            fprintf(stderr,
                "Particle %d outside of grid (%g, %g), marking as inactive.\n",
                i, job->particles[i].x, job->particles[i].y);
            job->particles[i].active = 0;
            continue;
        }

        /* Update particle element and mark element as occupied. */
        job->in_element[i] = p;
        job->elements[p].filled = 1;
        job->elements[p].n++;
        job->elements[p].m += job->particles[i].m;

    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void calculate_shapefunctions(job_t *job)
{
    int i, p, n;

    double xn;
    double yn;
    double xl;
    double yl;

    for (i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);
        p = job->in_element[i];
        n = job->elements[p].nodes[0];
        xn = job->nodes[n].x;
        yn = job->nodes[n].y;
        global_to_local_coords(&xl, &yl,
            job->particles[i].x, job->particles[i].y, 
            xn, yn, job->h);
        tent(&(job->h1[i]), &(job->h2[i]), &(job->h3[i]), &(job->h4[i]),
            xl, yl);
        grad_tent(
            &(job->b11[i]), &(job->b12[i]), &(job->b13[i]), &(job->b14[i]),
            &(job->b21[i]), &(job->b22[i]), &(job->b23[i]), &(job->b24[i]),
            xl, yl, job->h);
        if (xl < 0.0f || xl > 1.0f || yl < 0.0f || yl > 1.0f) {
            fprintf(stderr, "Particle %d outside of element %d (%g, %g).\n", i,
                p, xl, yl);
        }
        job->particles[i].xl = xl;
        job->particles[i].yl = yl;
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void calculate_node_velocity(job_t *job)
{
    int i;
    double m;

    for (i = 0; i < job->num_nodes; i++) {
        m = job->nodes[i].m;
        if (m > TOL) {
            job->nodes[i].x_t = job->nodes[i].mx_t / m;
            job->nodes[i].y_t = job->nodes[i].my_t / m;
        } else {
            job->nodes[i].x_t = 0;
            job->nodes[i].y_t = 0;
        }
    }

    for (i = 0; i < job->num_nodes; i++) {
        job->v_grid[NODAL_DOF * i + XDOF_IDX] = job->nodes[i].x_t;
        job->v_grid[NODAL_DOF * i + YDOF_IDX] = job->nodes[i].y_t;
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void calculate_strainrate(job_t *job)
{
    int i, j, k;
    int ce, nn[4];
/*    int method;*/
    double dx_tdy;
    double dy_tdx;
/*    char fname[] = "tmpx_sr.txt";*/

/*    FILE *fd;*/

/*    for (method = 0; method < 2; method++) { */ /* begin method loop */
/*    fname[3] = method + '0';*/
/*    fd = fopen(fname, "a");*/
/*    job->use_cpdi = method;*/
/*    fprintf(fd, "%g %d\n", job->t, job->use_cpdi);*/

    for (i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);
        job->particles[i].exx_t = 0;
        job->particles[i].exy_t = 0;
        job->particles[i].eyy_t = 0;
        job->particles[i].wxy_t = 0;

        if (job->use_cpdi) {
            /* loop over corners */
            for (j = 0; j < 4; j++) {
                ce = job->particles[i].corner_elements[j];
                /* corner is outside of particle domain. should probably deal with this better... */
                if (ce == -1) {
                    continue;
                }
                for (k = 0; k < 4; k++) {
                    nn[k] = job->elements[ce].nodes[k];

                    /* actual volume of particle here (not averaging volume) */
                    job->particles[i].exx_t += job->nodes[nn[k]].x_t * job->particles[i].grad_sc[j][k][S_XIDX];
                    job->particles[i].exy_t += 0.5 * (job->nodes[nn[k]].x_t * job->particles[i].grad_sc[j][k][S_YIDX] + job->nodes[nn[k]].y_t * job->particles[i].grad_sc[j][k][S_XIDX]);
                    job->particles[i].eyy_t += job->nodes[nn[k]].y_t * job->particles[i].grad_sc[j][k][S_YIDX];
                    job->particles[i].wxy_t += 0.5 * (job->nodes[nn[k]].x_t * job->particles[i].grad_sc[j][k][S_YIDX] - job->nodes[nn[k]].y_t * job->particles[i].grad_sc[j][k][S_XIDX]);
                }
            }
        } else {
            job->particles[i].exx_t = DX_N_TO_P(job, x_t, 1.0, i);
            job->particles[i].eyy_t = DY_N_TO_P(job, y_t, 1.0, i);

            dx_tdy = DY_N_TO_P(job, x_t, 1.0, i);
            dy_tdx = DX_N_TO_P(job, y_t, 1.0, i);

            job->particles[i].exy_t = 0.5 * (dx_tdy + dy_tdx);
            job->particles[i].wxy_t = 0.5 * (dx_tdy - dy_tdx);
        }
    }
/*    fprintf(fd, "%d %d: %g %g %g %g", method, i, job->particles[i].exx_t,*/
/*        job->particles[i].exy_t, job->particles[i].eyy_t, job->particles[i].wxy_t);*/

/*    fclose(fd);*/
    /*}*/ /*end method loop*/

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void update_stress(job_t *job)
{
    int i, j, k, method;
    int ce, nn[4];

    FILE *fd = fopen("tmp.txt", "a");

    for (method = 0; method < 2; method++) {

    job->use_cpdi = method;
    fprintf(fd, "%g %d\n", job->t, job->use_cpdi);

    for (i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);

        if (job->use_cpdi) {
            /* loop over corners */
            for (j = 0; j < 4; j++) {
                ce = job->particles[i].corner_elements[j];
                for (k = 0; k < 4; k++) {
                    nn[k] = job->elements[ce].nodes[k];
                    
                    /* actual volume of particle here (not averaging volume) */
                    job->nodes[nn[k]].fx += -job->particles[i].v * (
                    job->particles[i].grad_sc[j][k][S_XIDX] * job->particles[i].sxx + 
                    job->particles[i].grad_sc[j][k][S_YIDX] * job->particles[i].sxy);

                    job->nodes[nn[k]].fy += -job->particles[i].v * (
                    job->particles[i].grad_sc[j][k][S_XIDX] * job->particles[i].sxy + 
                    job->particles[i].grad_sc[j][k][S_YIDX] * job->particles[i].syy);
                }
            }
        } else {
            ACCUMULATE_WITH_MUL(fx, job, sxx, i, n, b1, -job->particles[i].v);
            ACCUMULATE_WITH_MUL(fx, job, sxy, i, n, b2, -job->particles[i].v);
            ACCUMULATE_WITH_MUL(fy, job, sxy, i, n, b1, -job->particles[i].v);
            ACCUMULATE_WITH_MUL(fy, job, syy, i, n, b2, -job->particles[i].v);
        }
    }

    for (i = 0; i < job->num_nodes; i++) {
        fprintf(fd, "f[%d] = [%g %g]^T\n", i, job->nodes[i].fx, job->nodes[i].fy);
    }

    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void update_internal_stress(job_t *job)
{
    int i, j, k, s;
    int p;
    int ce, nn[4];

    double sxx, sxy, syy;
    double pxx, pxy, pyx, pyy;
    double dudx, dudy, dvdx, dvdy, jdet;

    int method;

/*    char fname[] = "tmpx.txt";*/

/*    FILE *fd;*/

/*    for (method = 0; method < 2; method++) { */ /* begin method loop */
/*    fname[3] = method + '0';*/
/*    fd = fopen(fname, "a");*/
/*    job->use_cpdi = method;*/
/*    fprintf(fd, "%g %d\n", job->t, job->use_cpdi);*/

    for (i = 0; i < job->vec_len; i++) {
        job->f_int_grid[i] = 0;
    }

    for (i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);

        p = job->in_element[i];
        for (j = 0; j < 4; j++) {
            nn[j]  = job->elements[p].nodes[j];
        }

        if (job->use_cpdi) {
            /* loop over corners */
            for (j = 0; j < 4; j++) {
                ce = job->particles[i].corner_elements[j];
                /* corner is outside of particle domain. should probably deal with this better... */
                if (ce == -1) {
                    continue;
                }
                for (k = 0; k < 4; k++) {
                    nn[k] = job->elements[ce].nodes[k];

                    /* actual volume of particle here (not averaging volume) */
                    job->f_int_grid[NODAL_DOF * nn[k] + XDOF_IDX] += job->particles[i].v * (
                    job->particles[i].grad_sc[j][k][S_XIDX] * job->particles[i].sxx + 
                    job->particles[i].grad_sc[j][k][S_YIDX] * job->particles[i].sxy);

                    job->f_int_grid[NODAL_DOF * nn[k] + YDOF_IDX] += job->particles[i].v * (
                    job->particles[i].grad_sc[j][k][S_XIDX] * job->particles[i].sxy + 
                    job->particles[i].grad_sc[j][k][S_YIDX] * job->particles[i].syy);


/*                    printf("update_stress: %d\n", nn[k]);*/
/*                    printf("vol %g\ngrad_sc [%g %g]\ns [%g %g %g]\n", job->particles[i].v,*/
/*                     job->particles[i].grad_sc[j][k][S_XIDX], job->particles[i].grad_sc[j][k][S_YIDX], job->particles[i].sxx, job->particles[i].sxy, job->particles[i].syy);*/

/*                    for (s = 0; s < job->num_nodes; s++) {*/
/*                        printf("f[%d] = [%g %g]^T\n", s,*/
/*                            job->f_int_grid[NODAL_DOF * s + XDOF_IDX],*/
/*                            job->f_int_grid[NODAL_DOF * s + YDOF_IDX]);*/
/*                    }*/

                }
/*                exit(254);*/
            }
        } else {
#if 0
            dudx = job->u_grid[NODAL_DOF * nn[0] + XDOF_IDX]*job->b11[i];
            dudx += job->u_grid[NODAL_DOF * nn[1] + XDOF_IDX]*job->b12[i];
            dudx += job->u_grid[NODAL_DOF * nn[2] + XDOF_IDX]*job->b13[i];
            dudx += job->u_grid[NODAL_DOF * nn[3] + XDOF_IDX]*job->b14[i];

            dudy = job->u_grid[NODAL_DOF * nn[0] + XDOF_IDX]*job->b21[i];
            dudy += job->u_grid[NODAL_DOF * nn[1] + XDOF_IDX]*job->b22[i];
            dudy += job->u_grid[NODAL_DOF * nn[2] + XDOF_IDX]*job->b23[i];
            dudy += job->u_grid[NODAL_DOF * nn[3] + XDOF_IDX]*job->b24[i];

            dvdx = job->u_grid[NODAL_DOF * nn[0] + YDOF_IDX]*job->b11[i];
            dvdx += job->u_grid[NODAL_DOF * nn[1] + YDOF_IDX]*job->b12[i];
            dvdx += job->u_grid[NODAL_DOF * nn[2] + YDOF_IDX]*job->b13[i];
            dvdx += job->u_grid[NODAL_DOF * nn[3] + YDOF_IDX]*job->b14[i];

            dvdy = job->u_grid[NODAL_DOF * nn[0] + YDOF_IDX]*job->b21[i];
            dvdy += job->u_grid[NODAL_DOF * nn[1] + YDOF_IDX]*job->b22[i];
            dvdy += job->u_grid[NODAL_DOF * nn[2] + YDOF_IDX]*job->b23[i];
            dvdy += job->u_grid[NODAL_DOF * nn[3] + YDOF_IDX]*job->b24[i];

            jdet = 1 + dudx + dvdy;

            sxx = job->particles[i].sxx;
            sxy = job->particles[i].sxy;
            syy = job->particles[i].syy;

            pxx = jdet * ((1 - dudx) * sxx - dudy * sxy);
            pxy = jdet * ((-dvdx) * sxx + (1 - dvdy) * sxy);
            pyx = jdet * ((1 - dudx) * sxy - dudy * syy);
            pyy = jdet * ((-dvdx) * sxy + (1 - dvdy) * syy);
#endif

            sxx = job->particles[i].sxx;
            sxy = job->particles[i].sxy;
            syy = job->particles[i].syy;
            pxx = sxx;
            pxy = sxy;
            pyx = sxy;
            pyy = syy;

            job->f_int_grid[NODAL_DOF * nn[0] + XDOF_IDX] += job->particles[i].v*(job->b11[i]*pxx + job->b21[i]*pxy);
            job->f_int_grid[NODAL_DOF * nn[0] + YDOF_IDX] += job->particles[i].v*(job->b11[i]*pyx + job->b21[i]*pyy);
            job->f_int_grid[NODAL_DOF * nn[1] + XDOF_IDX] += job->particles[i].v*(job->b12[i]*pxx + job->b22[i]*pxy);
            job->f_int_grid[NODAL_DOF * nn[1] + YDOF_IDX] += job->particles[i].v*(job->b12[i]*pyx + job->b22[i]*pyy);
            job->f_int_grid[NODAL_DOF * nn[2] + XDOF_IDX] += job->particles[i].v*(job->b13[i]*pxx + job->b23[i]*pxy);
            job->f_int_grid[NODAL_DOF * nn[2] + YDOF_IDX] += job->particles[i].v*(job->b13[i]*pyx + job->b23[i]*pyy);
            job->f_int_grid[NODAL_DOF * nn[3] + XDOF_IDX] += job->particles[i].v*(job->b14[i]*pxx + job->b24[i]*pxy);
            job->f_int_grid[NODAL_DOF * nn[3] + YDOF_IDX] += job->particles[i].v*(job->b14[i]*pyx + job->b24[i]*pyy);
        }
    }

/*    for (i = 0; i < job->num_nodes; i++) {*/
/*        fprintf(fd, "f[%d] = [%g %g]^T\n", i,*/
/*            job->f_int_grid[NODAL_DOF * i + XDOF_IDX],*/
/*            job->f_int_grid[NODAL_DOF * i + YDOF_IDX]);*/
/*    }*/

/*    fclose(fd);*/
    /*}*/ /* end method loop*/


    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void implicit_solve(job_t *job)
{
    int i, j;
    int i_new, j_new;
    int r, s;
    int m, n;

    /* global and elemental indicies */
    int gi, gj;
    int ei, ej;
    double ke;

    /* lapack related*/
    int lda, ldb;

    /* dnrm2 */
    int incx = 1;
    double du_norm = 0;
    double du0_norm = 0;
    double q_norm = 0;
    double q0_norm = 0;

    /* For assembly of elements into global matrix. */
    int nn[4];

    /* iteration count */
    int k = 0;
    static int stable_timestep = 0;

    /* old timestep grid velocity and acceleration */
    double *v_grid_t;
    double *a_grid_t;

    /* mass weighted old timestep grid velocity and acceleration */
    double *mv_grid_t;
    double *ma_grid_t;
    double *m_grid_t;

    /* csparse */
    int nnz;
    int slda, sldb;
    cs *triplets;
    cs *smat;
    double *sb;
    int res = 0;

    /* for small timestep */
    double inv_dt = 1.0 / job->dt;
    double inv_dt_sq = inv_dt / job->dt;

    /* Timing code */
    struct timespec requestStart, requestEnd;
    long ns;

    lda = job->vec_len;
    ldb = lda;

    v_grid_t = (double *)malloc(lda * sizeof(double));
    a_grid_t = (double *)malloc(lda * sizeof(double));

    mv_grid_t = (double *)malloc(lda * sizeof(double));
    ma_grid_t = (double *)malloc(lda * sizeof(double));
    m_grid_t = (double *)malloc(lda * sizeof(double));

    /* keep old timestep velocity and acceleration */
    memcpy(v_grid_t, job->v_grid, lda * sizeof(double));
    memcpy(a_grid_t, job->a_grid, lda * sizeof(double));
/*    for (i = 0; i < ldb; i++) {*/
/*        mv_grid_t[i] = v_grid_t[i] * job->m_grid[i];*/
/*        ma_grid_t[i] = a_grid_t[i] * job->m_grid[i];*/
/*    }*/

    /* keep old stresses */
    for (i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);

        memcpy(job->particles[i].real_state, job->particles[i].state, DEPVAR * sizeof(double));
        job->particles[i].real_sxx = job->particles[i].sxx;
        job->particles[i].real_sxy = job->particles[i].sxy;
        job->particles[i].real_syy = job->particles[i].syy;
    }

start_implicit:
    k = 0;
    /* no initial grid deformation */
    for (i = 0; i < ldb; i++) {
        job->u_grid[i] = 0;
    }

    /* Generate and apply dirichlet boundary conditions given in BC file. */
    generate_dirichlet_bcs(job);

    /* Generate the new nodal numbers according to BCs. */
    generate_node_number_override(job);

    /*
        We have to change the velocities and accelerations if nodes are tied
        together. Weigh each component by nodal mass.
    */
    for (i = 0; i < lda; i++) {
        mv_grid_t[i] = 0;
        ma_grid_t[i] = 0;
        m_grid_t[i] = 0;
    }
    for (i = 0; i < lda; i++) {
        i_new = job->node_number_override[i];
        mv_grid_t[i_new] += v_grid_t[i] * job->m_grid[i];
        ma_grid_t[i_new] += a_grid_t[i] * job->m_grid[i];
        m_grid_t[i_new] += job->m_grid[i];
    }
    for (i = 0; i < lda; i++) {
        i_new = job->node_number_override[i];
        if (m_grid_t[i_new] == 0) {
            continue;
        }
        v_grid_t[i] = mv_grid_t[i_new] / m_grid_t[i_new];
        a_grid_t[i] = ma_grid_t[i_new] / m_grid_t[i_new];
    }

    /* Begin newton iterations */
    do {
        update_internal_stress(job);
        build_elemental_stiffness(job);

        /* starting timer */
        clock_gettime(CLOCK_REALTIME, &requestStart);

        /* build node and inverse maps, taking BCs into account. */
        for (i = 0; i < lda; i++) {
            job->node_u_map[i] = -1;
            job->inv_node_u_map[i] = -1;
        }
        j = 0;
        for (i = 0; i < lda; i++) {
            i_new = job->node_number_override[i];

            /* if the new node dof has already been assigned a place in the
                matrix, don't set it again */
            if (job->node_u_map[i_new] != -1) {
                continue;
            }

            /* otherwise, check if this node has any particles and the value
                has not been set by dirichlet BCs. */
            if (job->m_grid[i] > TOL && (job->u_dirichlet_mask[i] == 0)) {
                job->node_u_map[i_new] = j;
                j++;
            }
        }
        slda = j;   /* The number of free DOFs gives the sparse matrix size. */
        for (i = 0; i < lda; i++) {
            i_new = job->node_number_override[i];
            if ((j = job->node_u_map[i_new]) != -1) {
                job->inv_node_u_map[j] = i_new;
            }
        }

        /*  calculate number of nonzero elements (before summing duplicates). */
        nnz = 0;
        for (i = 0; i < job->num_elements; i++) {
            if (job->elements[i].filled != 0) {
                nnz += (NODAL_DOF * NODES_PER_ELEMENT) *
                        (NODAL_DOF * NODES_PER_ELEMENT);
                nnz += NODAL_DOF; /* for diagonal mass matrix entries. */
            }
        }

        /* allocate triplet array. */
        fprintf(stderr, "%d nonzeros (including duplicates) in K, which has size [%d x %d].\n", nnz, slda, slda);
        triplets = cs_spalloc(slda, slda, nnz, 1, 1);

        /* Calculate right hand side (load):
                Q = f_ext - f_int - M_g * (4 * u / dt^2 - 4 * v / dt - a) */
/*        memcpy(job->q_grid, job->m_grid, lda * sizeof(double));*/
/*        for (i = 0; i < lda; i++) {*/
/*            job->q_grid[i] = 0;*/
/*        }*/
/*        for (i = 0; i < lda; i++) {*/
/*            i_new = job->node_number_override[i];*/
/*            mv_grid_t[i_new] += v_grid_t[i] * job->m_grid[i];*/
/*            ma_grid_t[i_new] += a_grid_t[i] * job->m_grid[i];*/
/*            m_grid_t[i_new] += job->m_grid[i];*/
/*        }*/
/*        for (i = 0; i < lda; i++) {*/
/*            i_new = job->node_number_override[i];*/
/*            if (m_grid_t[i_new] == 0) {*/
/*                continue;*/
/*            }*/
/*            v_grid_t[i] = mv_grid_t[i_new] / m_grid_t[i_new];*/
/*            a_grid_t[i] = ma_grid_t[i_new] / m_grid_t[i_new];*/
/*        }*/
        for (i = 0; i < lda; i++) {
            i_new = job->node_number_override[i];
/*            job->q_grid[i_new] = (4 * mv_grid_t[i_new] * job->dt + ma_grid_t[i_new] * job->dt * job->dt);*/
            job->q_grid[i_new] = (1 * mv_grid_t[i_new] * job->dt);
        }
        for (i = 0; i < lda; i++) {
            i_new = job->node_number_override[i];
            /* XXX this is okay only because
                u_grid[i] == u_grid[i_new1] == u_grid[i_new2] etc ... */
/*            job->q_grid[i_new] += -job->m_grid[i] * (4 * job->u_grid[i]);*/
/*            job->q_grid[i_new] += (job->f_ext_grid[i] - job->f_int_grid[i]) * job->dt * job->dt;*/
            job->q_grid[i_new] += -job->m_grid[i] * (1 * job->u_grid[i]);
            job->q_grid[i_new] += (job->f_ext_grid[i] - job->f_int_grid[i]) * job->dt * job->dt;
        }

        /* assemble elements into global stiffness matrix */
        for (i = 0; i < job->num_elements; i++) {

            /* Ignore empty elements. */
            if (job->elements[i].filled == 0) {
                continue;
            }

            /* scale stiffness matrix */
            for (r = 0; r < (NODAL_DOF * NODES_PER_ELEMENT); r++) {
                for (s = 0; s < (NODAL_DOF * NODES_PER_ELEMENT); s++) {
                    job->elements[i].kku_element[r][s] *= (job->dt * job->dt);
                }
            }
            nn[0] = job->elements[i].nodes[0];
            nn[1] = job->elements[i].nodes[1];
            nn[2] = job->elements[i].nodes[2];
            nn[3] = job->elements[i].nodes[3];
            /*
                This presumes that the element level stiffness matrices are
                written as (x1, y1, c1, .... , x4, y4, c4) where the
                numbers indicate local node and the variable names are
                nodal degrees of freedom.
            */
            for (r = 0; r < NODES_PER_ELEMENT; r++) {
                for (m = 0; m < NODAL_DOF; m++) {
                    for (s = 0; s < NODES_PER_ELEMENT; s++) {
                        for (n = 0; n < NODAL_DOF; n++) {
                            gi = NODAL_DOF * nn[r] + m;
                            gj = NODAL_DOF * nn[s] + n;

                            /* renumber nodes as appropriate */
                            gi = job->node_number_override[gi];
                            gj = job->node_number_override[gj];

                            ei = NODAL_DOF * r + m;
                            ej = NODAL_DOF * s + n;
                            ke = job->elements[i].kku_element[ei][ej];

                            /*
                                If the global index is masked, don't add it
                                to the sparse matrix.
                            */
                            if (job->u_dirichlet_mask[gi] != 0) {
                                continue;
                            }

                            /* adjust load for dirichlet bcs. */
                            if (job->u_dirichlet_mask[gj] != 0) {
                                job->q_grid[gi] += (-job->u_dirichlet[gj] * ke);
                                continue;
                            }

                            res = cs_entry(triplets,
                                    job->node_u_map[gi],
                                    job->node_u_map[gj],
                                    ke);

                            if (res == 0) {
                                fprintf(stderr, "error adding stiffness entry\n");
                            }
                        }
                    }
                }
            }
        }

        for (i = 0; i < lda; i++) {
            i_new = job->node_number_override[i];
            if (job->node_u_map[i_new] != -1) {
                /* add diagonal mass matrix entries to
                    displacement dofs. */

/*                res = cs_entry(triplets,*/
/*                        job->node_u_map[i_new],*/
/*                        job->node_u_map[i_new],*/
/*                        4 * job->m_grid[i]);*/

                res = cs_entry(triplets,
                        job->node_u_map[i_new],
                        job->node_u_map[i_new],
                        1 * job->m_grid[i]);

                if (res == 0) {
                    fprintf(stderr, "error adding mass entry\n");
                }
            }
        }

        /* stopping timer */
        clock_gettime(CLOCK_REALTIME, &requestEnd);

        /* Calculate time it took */
    #define NS_PER_S 1E9
        ns = NS_PER_S * (requestEnd.tv_sec - requestStart.tv_sec )
          + ( requestEnd.tv_nsec - requestStart.tv_nsec );
        printf("Global Assembly[%d]: %lf s\n", k, (double)ns/NS_PER_S );
    #undef NS_PER_S

        /* reset stress state at beginning of newton iteration */
        for (i = 0; i < job->num_particles; i++) {
            CHECK_ACTIVE(job, i);
            
            memcpy(job->particles[i].state, job->particles[i].real_state, DEPVAR * sizeof(double));
            job->particles[i].sxx = job->particles[i].real_sxx;
            job->particles[i].sxy = job->particles[i].real_sxy;
            job->particles[i].syy = job->particles[i].real_syy;
        }

        sb = (double *)malloc(slda * sizeof(double));
        for (i = 0; i < slda; i++) {
            sb[i] = job->q_grid[job->inv_node_u_map[i]];
        }

        /* change to compressed format to sum up entries. */
        smat = cs_compress(triplets);

        /* sum entries with same indicies. */
        res = cs_dupl(smat);

        if (res == 0) {
            fprintf(stderr, "error summing entries\n");
            exit(255);
        }

        /* starting timer */
        clock_gettime(CLOCK_REALTIME, &requestStart);

        sldb = slda;
        if (k == 0) {
            q0_norm = dnrm2_(&sldb, sb, &incx);
        }
        q_norm = dnrm2_(&sldb, sb, &incx);

        if (!cs_lusol(1, smat, sb, 1e-12)) {
            fprintf(stderr, "lusol error!\n");
            if (cs_qrsol(1, smat, sb)) {
                fprintf(stderr, "qrsol error!\n");
                exit(255);
            }
        }

        if (k == 0) {
            du0_norm = dnrm2_(&sldb, sb, &incx);
        }
        du_norm = dnrm2_(&sldb, sb, &incx);

        /* stopping timer */
        clock_gettime(CLOCK_REALTIME, &requestEnd);

        /* Calculate time it took */
    #define NS_PER_S 1E9
        ns = NS_PER_S * (requestEnd.tv_sec - requestStart.tv_sec )
          + ( requestEnd.tv_nsec - requestStart.tv_nsec );
        printf("Iteration[%d]: Norm du: %f Norm q: %f Implicit Solve: %lf s\n", k, du_norm, q_norm, (double)ns/NS_PER_S );
    #undef NS_PER_S

        /* map solution back to entire grid */
        for (i = 0; i < ldb; i++) {
            i_new = job->node_number_override[i];
            if ((j = job->node_u_map[i_new]) != -1) {
                job->du_grid[i] = sb[j];
            } else {
                job->du_grid[i] = 0;
            }
        }

        /* apply dirichlet BCs to solution */
        for (i = 0; i < ldb; i++) {
            i_new = job->node_number_override[i];
            if (job->u_dirichlet_mask[i_new] != 0) {
                job->du_grid[i] = job->u_dirichlet[i_new];
            }
        }

        for (i = 0; i < ldb; i++) {
            /* update grid displacement */
            job->u_grid[i] += job->du_grid[i];

            /* update grid velocity */
/*            job->v_grid[i] = 2 * (job->u_grid[i] / job->dt) - v_grid_t[i];*/
/*            job->v_grid[i] = (job->u_grid[i] / job->dt) - v_grid_t[i];*/
            job->v_grid[i] = (job->u_grid[i] / job->dt);
        }

        for (i = 0; i < job->num_nodes; i++) {
            job->nodes[i].x_t = job->v_grid[NODAL_DOF * i + XDOF_IDX];
            job->nodes[i].y_t = job->v_grid[NODAL_DOF * i + YDOF_IDX];
        }

        calculate_strainrate(job);
        calculate_stress(job);

        cs_spfree(triplets);
        cs_spfree(smat);
        free(sb);
        k++;

        if (du_norm < job->implicit.du_norm_converged) {
            printf("norm of du  = %e: Accepting as converged.\n", du_norm);
            break;
        }

        if (q_norm*du_norm > 10*(q0_norm*du0_norm)
                || du_norm > 10*du0_norm
                || k > job->implicit.unstable_iteration_count) {
            job->dt = job->dt * 0.8;
            stable_timestep = 0;
            printf("Trouble converging, modifying Timestep to %f.\n", job->dt);
            goto start_implicit;
        }

        if (job->dt < job->timestep.dt_min) {
            fprintf(stderr, "Timestep %g is too small, aborting.\n", job->dt);
            exit(-1);
        }

    } while ((du_norm / du0_norm) > job->implicit.du_norm_ratio ||
        ((du_norm*q_norm) / (du0_norm*q0_norm)) > job->implicit.q_norm_ratio);

    stable_timestep++;
    if (job->timestep.allow_dt_increase != 0
        && stable_timestep > job->timestep.stable_dt_threshold 
        && job->dt < job->timestep.dt_max) {
        job->dt = job->dt * 1.25;
        if (job->dt > job->timestep.dt_max) {
            job->dt = job->timestep.dt_max;
        }
        printf("Increasing timestep to %f.\n", job->dt);
        stable_timestep = 0;
    }

    /* update grid acceleration */
    for (i = 0; i < ldb; i++) {
/*        job->a_grid[i] = (4 * job->u_grid[i] * inv_dt_sq - 4 * v_grid_t[i] * inv_dt - a_grid_t[i]);*/
        job->a_grid[i] = (job->u_grid[i] * inv_dt_sq - v_grid_t[i] * inv_dt);
    }

    /* repackage for use with macros */
    for (i = 0; i < job->num_nodes; i++) {
        job->nodes[i].x_tt = job->a_grid[NODAL_DOF * i + XDOF_IDX];
        job->nodes[i].y_tt = job->a_grid[NODAL_DOF * i + YDOF_IDX];

        job->nodes[i].x_t = job->v_grid[NODAL_DOF * i + XDOF_IDX];
        job->nodes[i].y_t = job->v_grid[NODAL_DOF * i + YDOF_IDX];

        job->nodes[i].ux = job->u_grid[NODAL_DOF * i + XDOF_IDX];
        job->nodes[i].uy = job->u_grid[NODAL_DOF * i + YDOF_IDX];
    }

    free(v_grid_t);
    free(a_grid_t);
    free(mv_grid_t);
    free(ma_grid_t);
    free(m_grid_t);

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void map_to_grid(job_t *job)
{
    int i;

    for (i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);

        /* Mass. */
        ACCUMULATE(m,job,m,i,n,h);

        /* Momentum. */
        ACCUMULATE_WITH_MUL(mx_t,job,x_t,i,n,h,job->particles[i].m);
        ACCUMULATE_WITH_MUL(my_t,job,y_t,i,n,h,job->particles[i].m);

        /* Pseudoforce. */
        ACCUMULATE_WITH_MUL(mx_tt,job,x_tt,i,n,h,job->particles[i].m);
        ACCUMULATE_WITH_MUL(my_tt,job,y_tt,i,n,h,job->particles[i].m);

        /* Body forces. */
        ACCUMULATE_WITH_MUL(fx,job,bx,i,n,h,job->particles[i].m);
        ACCUMULATE_WITH_MUL(fy,job,by,i,n,h,job->particles[i].m);
    }

    /* for implicit method, map to another array */
    for (i = 0; i < job->num_nodes; i++) {
        /* used lumped mass */
        job->m_grid[NODAL_DOF * i + XDOF_IDX] = job->nodes[i].m;
        job->m_grid[NODAL_DOF * i + YDOF_IDX] = job->nodes[i].m;

        /* need previous timestep's acceleration for implicit method */
        if (job->nodes[i].m > TOL) {
            job->a_grid[NODAL_DOF * i + XDOF_IDX] = job->nodes[i].mx_tt / job->nodes[i].m;
            job->a_grid[NODAL_DOF * i + YDOF_IDX] = job->nodes[i].my_tt  / job->nodes[i].m;
        } else {
            job->a_grid[NODAL_DOF * i + XDOF_IDX] = 0;
            job->a_grid[NODAL_DOF * i + YDOF_IDX] = 0;
        }

        /* external forces */
        job->f_ext_grid[NODAL_DOF * i + XDOF_IDX] = job->nodes[i].fx;
        job->f_ext_grid[NODAL_DOF * i + YDOF_IDX] = job->nodes[i].fy;
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void move_particles(job_t *job)
{
    int i;
    double dux, a_x_t;
    double duy, a_y_t;

    for (i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);

        dux = (N_TO_P(job, ux, i));
        duy = (N_TO_P(job, uy, i));

        job->particles[i].x += dux;
        job->particles[i].y += duy;
        job->particles[i].ux += dux;
        job->particles[i].uy += duy;

        while(job->particles[i].x < 0) { job->particles[i].x += 1.0; }
        while(job->particles[i].x > 1) { job->particles[i].x -= 1.0; }

        a_x_t = job->particles[i].x_tt;
        a_y_t = job->particles[i].y_tt;

        job->particles[i].x_tt = (N_TO_P(job, x_tt, i));
        job->particles[i].y_tt = (N_TO_P(job, y_tt, i));

/*        job->particles[i].x_t += 0.5 * job->dt * (job->particles[i].x_tt + a_x_t);*/
/*        job->particles[i].y_t += 0.5 * job->dt * (job->particles[i].y_tt + a_y_t);*/
        job->particles[i].x_t = (N_TO_P(job, x_t, i));
        job->particles[i].y_t = (N_TO_P(job, y_t, i));
    }
    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void update_deformation_gradient(job_t *job)
{
    int i;

    double delta;
    double delta_inv;
    double m11, m12, m21, m22;
    double a, b, c, d;
    double cd, sd, f;
    double dsq;
    double tmp_Fxx;
    double tmp_Fxy;
    double tmp_Fyx;
    double tmp_Fyy;

    for (i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);

        a = job->dt * DX_N_TO_P(job, x_t, 1, i);
        b = job->dt * DY_N_TO_P(job, x_t, 1, i);
        c = job->dt * DX_N_TO_P(job, y_t, 1, i);
        d = job->dt * DY_N_TO_P(job, y_t, 1, i);

        dsq = (a - d) * (a - d) + 4 * b * c;
        delta = 0.5 * sqrt(abs(dsq));
        delta_inv = (1.0f / delta);
        if (dsq == 0 || delta < 1e-15) {
            cd = 1;
            sd = 1;
            delta_inv = 1;
        } else if (dsq > 0) {
            cd = cosh(delta);
            sd = sinh(delta);
        } else {
            cd = cos(delta);
            sd = sin(delta);
        }

        m11 = cd + 0.5 * (a - d) * sd * delta_inv;
        m12 = b * sd * delta_inv;
        m21 = c * sd * delta_inv;
        m22 = cd - 0.5 * (a - d) * sd * delta_inv;

        f = exp(0.5 * (a + d));
        tmp_Fxx = f * (m11);
        tmp_Fxy = f * (m12);
        tmp_Fyx = f * (m21);
        tmp_Fyy = f * (m22);

        job->particles[i].Fxx = job->particles[i].Fxx * tmp_Fxx + job->particles[i].Fyx * tmp_Fxy;
        job->particles[i].Fxy = job->particles[i].Fxy * tmp_Fxx + job->particles[i].Fyy * tmp_Fxy;
        job->particles[i].Fyx = job->particles[i].Fxx * tmp_Fyx + job->particles[i].Fyx * tmp_Fyy;
        job->particles[i].Fyy = job->particles[i].Fxy * tmp_Fyx + job->particles[i].Fyy * tmp_Fyy;

/*            job->particles[i].v = (job->particles[i].Fxx * job->particles[i].Fyy - job->particles[i].Fxy * job->particles[i].Fyx) * job->particles[i].v0;*/
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void update_particle_vectors(job_t *job)
{
    int i;
    for (i = 0; i < job->num_particles; i++) {
        job->particles[i].r1[0] =
            job->particles[i].Fxx * job->particles[i].r1_initial[0] +
            job->particles[i].Fxy * job->particles[i].r1_initial[1];
        job->particles[i].r1[1] =
            job->particles[i].Fyx * job->particles[i].r1_initial[0] +
            job->particles[i].Fyy * job->particles[i].r1_initial[1];
        job->particles[i].r2[0] =
            job->particles[i].Fxx * job->particles[i].r2_initial[0] +
            job->particles[i].Fxy * job->particles[i].r2_initial[1];
        job->particles[i].r2[1] =
            job->particles[i].Fyx * job->particles[i].r2_initial[0] +
            job->particles[i].Fyy * job->particles[i].r2_initial[1];
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void update_corner_positions(job_t *job)
{
    int i, j, k;
    int e, n;

/*    FILE *fd;*/
/*    fd = fopen("pnodes.txt", "a");*/

/*    for (i = 0; i < (job->num_particles * job->num_nodes); i++) {*/
/*        job->phi[i] = 0;*/
/*    }*/

    for (i = 0; i < job->num_particles; i++) {
        /* Find corner positions from (adjusted) vectors. */
        job->particles[i].corners[2][S_XIDX] = job->particles[i].x +
            job->particles[i].r1[S_XIDX] + job->particles[i].r2[S_XIDX];
        job->particles[i].corners[2][S_YIDX] = job->particles[i].y +
            job->particles[i].r1[S_YIDX] + job->particles[i].r2[S_YIDX];

        job->particles[i].corners[3][S_XIDX] = job->particles[i].x -
            job->particles[i].r1[S_XIDX] + job->particles[i].r2[S_XIDX];
        job->particles[i].corners[3][S_YIDX] = job->particles[i].y -
            job->particles[i].r1[S_YIDX] + job->particles[i].r2[S_YIDX];

        job->particles[i].corners[0][S_XIDX] = job->particles[i].x -
            job->particles[i].r1[S_XIDX] - job->particles[i].r2[S_XIDX];
        job->particles[i].corners[0][S_YIDX] = job->particles[i].y -
            job->particles[i].r1[S_YIDX] - job->particles[i].r2[S_YIDX];

        job->particles[i].corners[1][S_XIDX] = job->particles[i].x +
            job->particles[i].r1[S_XIDX] - job->particles[i].r2[S_XIDX];
        job->particles[i].corners[1][S_YIDX] = job->particles[i].y +
            job->particles[i].r1[S_YIDX] - job->particles[i].r2[S_YIDX];

        for (j = 0; j < 4; j++) {
            /* Figure out which element each corner is in. */
            job->particles[i].corner_elements[j] = WHICH_ELEMENT(job->particles[i].corners[j][S_XIDX], job->particles[i].corners[j][S_YIDX], job->N, job->h);

            e = job->particles[i].corner_elements[j];
            n = job->elements[e].nodes[0];

/*            fprintf(fd, "%d:", i);*/
/*            for (k = 0; k < 4; k++) {*/
/*                fprintf(fd, " %d", job->elements[e].nodes[k]);*/
/*            }*/
/*            fprintf(fd, "\n");*/

            /* calculate local coordinates of corners */
            if (job->particles[i].corner_elements[j] != -1) {
                global_to_local_coords(
                    &(job->particles[i].cornersl[j][S_XIDX]),
                    &(job->particles[i].cornersl[j][S_YIDX]),
                    job->particles[i].corners[j][S_XIDX],
                    job->particles[i].corners[j][S_YIDX],
                    job->nodes[n].x,
                    job->nodes[n].y,
                    job->h
                );

                tent(
                    &(job->particles[i].sc[j][0]),
                    &(job->particles[i].sc[j][1]),
                    &(job->particles[i].sc[j][2]),
                    &(job->particles[i].sc[j][3]),
                    job->particles[i].cornersl[j][S_XIDX],
                    job->particles[i].cornersl[j][S_YIDX]
                );
            }
        }

        for (j = 0; j < 4; j++) {
            for (k = 0; k < 4; k++) {

                switch(k) {
                    case 0:
                        job->particles[i].grad_sc[k][j][S_XIDX] = (job->particles[i].sc[0][j])*(job->particles[i].r1[S_YIDX] - job->particles[i].r2[S_YIDX]);

                        job->particles[i].grad_sc[k][j][S_YIDX] = (job->particles[i].sc[0][j])*(job->particles[i].r2[S_XIDX] - job->particles[i].r1[S_XIDX]);
                    break;
                    case 1:
                        job->particles[i].grad_sc[k][j][S_XIDX] = (job->particles[i].sc[1][j])*(job->particles[i].r1[S_YIDX] + job->particles[i].r2[S_YIDX]);

                        job->particles[i].grad_sc[k][j][S_YIDX] = (job->particles[i].sc[1][j])*(-job->particles[i].r1[S_XIDX] - job->particles[i].r2[S_XIDX]);
                    break;
                    case 2:
                        job->particles[i].grad_sc[k][j][S_XIDX] = (-job->particles[i].sc[2][j])*(job->particles[i].r1[S_YIDX] - job->particles[i].r2[S_YIDX]);

                        job->particles[i].grad_sc[k][j][S_YIDX] = (-job->particles[i].sc[2][j])*(job->particles[i].r2[S_XIDX] - job->particles[i].r1[S_XIDX]);
                    break;
                    case 3:
                        job->particles[i].grad_sc[k][j][S_XIDX] = (-job->particles[i].sc[3][j])*(job->particles[i].r1[S_YIDX] + job->particles[i].r2[S_YIDX]);

                        job->particles[i].grad_sc[k][j][S_YIDX] = (-job->particles[i].sc[3][j])*(-job->particles[i].r1[S_XIDX] - job->particles[i].r2[S_XIDX]);
                    break;
                }
/*                job->particles[i].grad_sc[k][j][S_XIDX] = (job->particles[i].sc[0][j] - job->particles[i].sc[2][j])*(job->particles[i].r1[S_YIDX] - job->particles[i].r2[S_YIDX]) +*/
/*                (job->particles[i].sc[1][j] - job->particles[i].sc[3][j])*(job->particles[i].r1[S_YIDX] + job->particles[i].r2[S_YIDX]);*/

/*                job->particles[i].grad_sc[k][j][S_YIDX] = (job->particles[i].sc[0][j] - job->particles[i].sc[2][j])*(job->particles[i].r2[S_XIDX] - job->particles[i].r1[S_XIDX]) +*/
/*                (job->particles[i].sc[1][j] - job->particles[i].sc[3][j])*(-job->particles[i].r1[S_XIDX] - job->particles[i].r2[S_XIDX]);*/

/*                printf("%d,%d sc [%g %g %g %g]\n", j, k, job->particles[i].sc[k][0], job->particles[i].sc[k][1], job->particles[i].sc[k][2], job->particles[i].sc[k][3]);*/
/*                printf("%d,%d sc [%g %g]\n", j, k, job->particles[i].grad_sc[k][j][S_XIDX],*/
/*                    job->particles[i].grad_sc[k][j][S_YIDX]);*/

                /* should be averaging domain, but use particle volume for now */
                job->particles[i].grad_sc[k][j][S_XIDX] *= (1.0 / job->particles[i].v);
                job->particles[i].grad_sc[k][j][S_YIDX] *= (1.0 / job->particles[i].v);

                if (job->particles[i].corner_elements[k] == -1) {
                    continue;
                }

                /* phi rows are nodal indicies, cols are particle indicies */
/*                job->phi[job->num_particles * job->elements[job->particles[i].corner_elements[k]].nodes[j] + i] += 0.25 * job->particles[i].sc[k][j];*/
            }
        }

    }

/*    for (i = 0; i < job->num_nodes; i++) {*/
/*        for (j = 0; j < job->num_particles; j++) {*/
/*            fprintf(fd, " %g", job->phi[job->num_particles * i + j]);*/
/*        }*/
/*        fprintf(fd, "\n");*/
/*    }*/
/*    fprintf(fd, "\n");*/

/*    fclose(fd);*/

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void update_particle_densities(job_t *job)
{
    int i;
    int j;

    double left, right, top, bottom;
    double h_spacing, v_spacing;

    for (i = 0; i < job->num_elements; i++) {
        job->elements[i].filled = 0;
        job->elements[i].n = 0;
        job->elements[i].m = 0;
    }

    for (i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);
        job->elements[job->in_element[i]].filled = 1;
        job->elements[job->in_element[i]].n++;
        job->elements[job->in_element[i]].m += job->particles[i].m;
        
    }

    for (i = 0; i < job->num_nodes; i++) {
        job->nodes[i].num_filled_element_neighbors = 0;
        job->nodes[i].mass_filled_element_neighbors = 0;
    }

    for (i = 0; i < job->num_elements; i++) {
        if (job->elements[i].filled) {
            for (j = 0; j < 4; j++) {
                job->nodes[job->elements[i].nodes[j]].num_filled_element_neighbors++;
                job->nodes[job->elements[i].nodes[j]].mass_filled_element_neighbors += job->elements[i].m;
            }
        }
    }

    for (i = 0; i < job->num_elements; i++) {
        h_spacing = 2.0f;
        v_spacing = 2.0f;

        if (job->elements[i].neighbors[0] == -1) {
            right = job->elements[i].n;
            h_spacing = 1.0f;
        } else {
            right = job->elements[job->elements[i].neighbors[0]].n;
        }

        if (job->elements[i].neighbors[4] == -1) {
            left = job->elements[i].n;
            h_spacing = 1.0f;
        } else {
            left = job->elements[job->elements[i].neighbors[4]].n;
        }

        if (job->elements[i].neighbors[2] == -1) {
            top = job->elements[i].n;
            v_spacing = 1.0f;
        } else {
            top = job->elements[job->elements[i].neighbors[2]].n;
        }

        if (job->elements[i].neighbors[6] == -1) {
            bottom = job->elements[i].n;
            v_spacing = 1.0f;
        } else {
            bottom = job->elements[job->elements[i].neighbors[6]].n;
        }

        job->elements[i].grad_x = (right - left) / h_spacing;
        job->elements[i].grad_y = (top - bottom) / v_spacing;

        job->elements[i].grad_mag = hypot(job->elements[i].grad_x, job->elements[i].grad_y);

        /* point normal outward from body (negative of gradient). */
        if (job->elements[i].grad_mag > 0) {
            job->elements[i].n_x = -job->elements[i].grad_x / job->elements[i].grad_mag;
            job->elements[i].n_y = -job->elements[i].grad_y / job->elements[i].grad_mag;
            job->elements[i].n_theta = atan2(job->elements[i].n_y, job->elements[i].n_x);
        }

        /* printing */
/*        if (job->elements[i].grad_mag > 1.0f) {*/
/*            printf("%d: %f %f %f %f\n", i, left, right, top, bottom);*/
/*            printf("%d: %f\n", i, job->elements[i].grad_mag);*/
/*            printf("%d: %f %f\n", i, job->elements[i].n_x, job->elements[i].n_y);*/
/*        }*/

    }

        #define GRAD_THRESHOLD 0.50f

    for (i = 0; i < job->num_nodes; i++) {
        job->nodes[i].rho = job->nodes[i].mass_filled_element_neighbors / (job->nodes[i].num_filled_element_neighbors * job->h * job->h);
    }

    for (i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);

        /* Fixed version */
/*        if (job->elements[job->in_element[i]].grad_mag > GRAD_THRESHOLD) {*/
/*            job->particles[i].v = job->particles[i].v * exp (job->dt * (job->particles[i].exx_t + job->particles[i].eyy_t));*/
/*        } else {*/
/*            job->particles[i].v = job->h * job->h / (job->elements[job->in_element[i]].n);*/
/*        }*/

/*        if (job->use_cpdi) {*/
            job->particles[i].v = job->particles[i].v0 *
                ((job->particles[i].Fxx * job->particles[i].Fyy) - 
                job->particles[i].Fxy * job->particles[i].Fyx);
/*        } else {*/
/*            job->particles[i].v = job->h * job->h / (job->elements[job->in_element[i]].n);*/
/*        }*/
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void mpm_cleanup(job_t *job)
{
    free(job->nodes);
    free(job->particles);
    free(job->elements);
    
    free(job->in_element);

    free(job->h1);
    free(job->h2);
    free(job->h3);
    free(job->h4);
/*    free(job->h5);*/
/*    free(job->h6);*/
/*    free(job->h7);*/
/*    free(job->h8);*/
/*    free(job->h9);*/

    free(job->b11);
    free(job->b12);
    free(job->b13);
    free(job->b14);
/*    free(job->b15);*/
/*    free(job->b16);*/
/*    free(job->b17);*/
/*    free(job->b18);*/
/*    free(job->b19);*/

    free(job->b21);
    free(job->b22);
    free(job->b23);
    free(job->b24);
/*    free(job->b25);*/
/*    free(job->b26);*/
/*    free(job->b27);*/
/*    free(job->b28);*/
/*    free(job->b29);*/

    free(job);

    return;
}
/*----------------------------------------------------------------------------*/

