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
#include "particle.h"
#include "node.h"
#include "process.h"
#include "material.h"
#include <suitesparse/cs.h>

/*#include <omp.h>*/
#include <signal.h>

#define TOL 1e-10

#define signum(x) ((int)((0 < x) - (x < 0)))

#define ijton(i,j,N) ((j)*N + (i))
#define linidx_fortran(i,j,C) ((i)*C + (j))
#define N_TO_P(j,tok,i) SMEAR(j,tok,i,h)
#define DX_N_TO_P(j,tok,c,i) (c) * SMEAR(j,tok,i,b1)
#define DY_N_TO_P(j,tok,c,i) (c) * SMEAR(j,tok,i,b2)

#define IS_VALID_ELEMENT_COORD4(r,c,N) \
    (((r) < (N)) && ((r) >= 0) && ((c) < (N)) && ((c) >= 0))

#define FILL_ELEMENT_NEIGHBOR(en,r,c,N) \
    if (IS_VALID_ELEMENT_COORD((r),(c),(N))) { \
        en = (r)*(N) + (c); \
    } else { \
        en = -1; \
    }

/*#define QUAD_ELEMENTS*/

/* trololololol */
#define h1(j,i) j->particles[i].h[0]
#define h2(j,i) j->particles[i].h[1]
#define h3(j,i) j->particles[i].h[2]
#define h4(j,i) j->particles[i].h[3]
#define h5(j,i) j->particles[i].h[4]
#define h6(j,i) j->particles[i].h[5]
#define h7(j,i) j->particles[i].h[6]
#define h8(j,i) j->particles[i].h[7]
#define h9(j,i) j->particles[i].h[8]

#ifdef QUAD_ELEMENTS
    #define SMEAR SMEAR9
    #define ACCUMULATE ACCUMULATE9
    #define ACCUMULATE_WITH_MUL ACCUMULATE_WITH_MUL9
    #define WHICH_ELEMENT WHICH_ELEMENT9
#else
    #define SMEAR SMEAR4
    #define ACCUMULATE ACCUMULATE4
    #define ACCUMULATE_WITH_MUL ACCUMULATE_WITH_MUL4
    #define WHICH_ELEMENT WHICH_ELEMENT4
    #define IS_VALID_ELEMENT_COORD IS_VALID_ELEMENT_COORD4
#endif

#define __E(j,p) j->in_element[p]
#define __N(j,p,n) j->elements[__E(j,p)].nodes[n]
#define __NE(j,e,n) j->elements[e].nodes[n]

#define WHICH_ELEMENT4(xp,yp,N,h) \
    (floor(xp/h) + floor(yp/h)*(N-1))

#define WHICH_ELEMENT9(xp,yp,N,h) \
    (floor(xp/(2*h)) + floor(yp/(2*h))*(N-1)/2)

#define ACCUMULATE4(acc_tok,j,tok,i,n,s) \
    j->nodes[__N(j,i,0)].acc_tok += j->s ## 1[i] * j->particles[i].tok; \
    j->nodes[__N(j,i,1)].acc_tok += j->s ## 2[i] * j->particles[i].tok; \
    j->nodes[__N(j,i,2)].acc_tok += j->s ## 3[i] * j->particles[i].tok; \
    j->nodes[__N(j,i,3)].acc_tok += j->s ## 4[i] * j->particles[i].tok;

#define ACCUMULATE9(acc_tok,j,tok,i,n,s) \
    ACCUMULATE4(acc_tok,j,tok,i,n,s); \
    j->nodes[__N(j,i,4)].acc_tok += j->s ## 5[i] * j->particles[i].tok; \
    j->nodes[__N(j,i,5)].acc_tok += j->s ## 6[i] * j->particles[i].tok; \
    j->nodes[__N(j,i,6)].acc_tok += j->s ## 7[i] * j->particles[i].tok; \
    j->nodes[__N(j,i,7)].acc_tok += j->s ## 8[i] * j->particles[i].tok; \
    j->nodes[__N(j,i,8)].acc_tok += j->s ## 9[i] * j->particles[i].tok;

#define ACCUMULATE_WITH_MUL4(acc_tok,j,tok,i,n,s,c) \
    j->nodes[__N(j,i,0)].acc_tok += j->s ## 1[i] * j->particles[i].tok * (c); \
    j->nodes[__N(j,i,1)].acc_tok += j->s ## 2[i] * j->particles[i].tok * (c); \
    j->nodes[__N(j,i,2)].acc_tok += j->s ## 3[i] * j->particles[i].tok * (c); \
    j->nodes[__N(j,i,3)].acc_tok += j->s ## 4[i] * j->particles[i].tok * (c);

#define ACCUMULATE_WITH_MUL9(acc_tok,j,tok,i,n,s,c) \
    ACCUMULATE_WITH_MUL4(acc_tok,j,tok,i,n,s,c); \
    j->nodes[__N(j,i,4)].acc_tok += j->s ## 5[i] * j->particles[i].tok * (c); \
    j->nodes[__N(j,i,5)].acc_tok += j->s ## 6[i] * j->particles[i].tok * (c); \
    j->nodes[__N(j,i,6)].acc_tok += j->s ## 7[i] * j->particles[i].tok * (c); \
    j->nodes[__N(j,i,7)].acc_tok += j->s ## 8[i] * j->particles[i].tok * (c); \
    j->nodes[__N(j,i,8)].acc_tok += j->s ## 9[i] * j->particles[i].tok * (c);

#define SMEAR4(j,tok,i,s) ( \
    j->s ## 1[i] * j->nodes[__N(j,i,0)].tok + \
    j->s ## 2[i] * j->nodes[__N(j,i,1)].tok + \
    j->s ## 3[i] * j->nodes[__N(j,i,2)].tok + \
    j->s ## 4[i] * j->nodes[__N(j,i,3)].tok \
)

#define SMEAR9(j,tok,i,s) ( \
    SMEAR4(j,tok,i,s) + \
    j->s ## 5[i] * j->nodes[__N(j,i,4)].tok + \
    j->s ## 6[i] * j->nodes[__N(j,i,5)].tok + \
    j->s ## 7[i] * j->nodes[__N(j,i,6)].tok + \
    j->s ## 8[i] * j->nodes[__N(j,i,7)].tok + \
    j->s ## 9[i] * j->nodes[__N(j,i,8)].tok \
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
    job->particles = (particle_t *)malloc(num_particles * sizeof(particle_t));
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

    /* Get node coordinates. */
    job->nodes = (node_t *)malloc(job->num_nodes * sizeof(node_t));
    for (i = 0; i < job->num_nodes; i++) {
        node_number_to_coords(&(job->nodes[i].x), &(job->nodes[i].y), i, N, h);
        /* printf("preprocessing: xn=%g yn=%g\n", job->nodes[i].x, job->nodes[i].y); */
    }

    /* Get node numbering for elements. */
    job->elements = (element_t *)malloc(job->num_elements * sizeof(element_t));
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
    job->h5 = (double *)malloc(job->num_particles * sizeof(double));
    job->h6 = (double *)malloc(job->num_particles * sizeof(double));
    job->h7 = (double *)malloc(job->num_particles * sizeof(double));
    job->h8 = (double *)malloc(job->num_particles * sizeof(double));
    job->h9 = (double *)malloc(job->num_particles * sizeof(double));

    job->b11 = (double *)malloc(job->num_particles * sizeof(double));
    job->b12 = (double *)malloc(job->num_particles * sizeof(double));
    job->b13 = (double *)malloc(job->num_particles * sizeof(double));
    job->b14 = (double *)malloc(job->num_particles * sizeof(double));
    job->b15 = (double *)malloc(job->num_particles * sizeof(double));
    job->b16 = (double *)malloc(job->num_particles * sizeof(double));
    job->b17 = (double *)malloc(job->num_particles * sizeof(double));
    job->b18 = (double *)malloc(job->num_particles * sizeof(double));
    job->b19 = (double *)malloc(job->num_particles * sizeof(double));

    job->b21 = (double *)malloc(job->num_particles * sizeof(double));
    job->b22 = (double *)malloc(job->num_particles * sizeof(double));
    job->b23 = (double *)malloc(job->num_particles * sizeof(double));
    job->b24 = (double *)malloc(job->num_particles * sizeof(double));
    job->b25 = (double *)malloc(job->num_particles * sizeof(double));
    job->b26 = (double *)malloc(job->num_particles * sizeof(double));
    job->b27 = (double *)malloc(job->num_particles * sizeof(double));
    job->b28 = (double *)malloc(job->num_particles * sizeof(double));
    job->b29 = (double *)malloc(job->num_particles * sizeof(double));

    #define DIMENSION 2
    /* max size of u_grid is dim * N * N. */
    job->vec_len = DIMENSION * job->N * job->N;
    job->u_grid = (double *)malloc(DIMENSION * sizeof(double) * job->N * job->N);
    job->du_grid = (double *)malloc(DIMENSION * sizeof(double) * job->N * job->N);
    job->v_grid = (double *)malloc(DIMENSION * sizeof(double) * job->N * job->N);
    job->a_grid = (double *)malloc(DIMENSION * sizeof(double) * job->N * job->N);
    job->m_grid = (double *)malloc(DIMENSION * sizeof(double) * job->N * job->N);
    job->f_ext_grid = (double *)malloc(DIMENSION * sizeof(double) * job->N * job->N);
    job->f_int_grid = (double *)malloc(DIMENSION * sizeof(double) * job->N * job->N);
    job->q_grid = (double *)malloc(DIMENSION * sizeof(double) * job->N * job->N);
    job->node_u_map = (int *)malloc(DIMENSION * sizeof(int) * job->N * job->N);
    job->inv_node_u_map = (int *)malloc(DIMENSION * sizeof(int) * job->N * job->N);

    job->kku_grid = (double *)malloc(
                        (DIMENSION * job->N * job->N) *
                        (DIMENSION * sizeof(double) * job->N * job->N));
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
    for (i = 0; i < (job->vec_len * job->vec_len); i++) {
        job->kku_grid[i] = 0;
    }
    #undef DIMENSION

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
    }

    /* intialize state variables */
    material_init(job);

    /* Set default timestep. */
    job->dt = 0.4 * job->h * sqrt(job->particles[0].m/(job->particles[0].v * EMOD));

    /* By default, don't output energy data. */
    job->ke_data = NULL;

    return job;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void mpm_step(job_t *job)
{
    int i;
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

        /* first component is x, second is y */
        job->u_grid[2 * i] = 0;
        job->u_grid[2 * i + 1] = 0;

        job->du_grid[2 * i] = 0;
        job->du_grid[2 * i + 1] = 0;

        job->f_int_grid[2 * i] = 0;
        job->f_int_grid[2 * i + 1] = 0;

        job->f_ext_grid[2 * i] = 0;
        job->f_ext_grid[2 * i + 1] = 0;
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

    implicit_solve(job);

    /* Update particle position and velocity. */
    move_particles(job);

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

        /* Update particle element and mark element as occupied. */
        job->in_element[i] = p;
        job->elements[job->in_element[i]].filled = 1;
        job->elements[job->in_element[i]].n++;
        job->elements[job->in_element[i]].m += job->particles[i].m;

        /* XXX ugly hack for domain size... */
        if (job->particles[i].x < 0.0 || job->particles[i].x > 1.0
            || job->particles[i].y < 0.0 || job->particles[i].y > 1.0) {
            fprintf(stderr,
                "Particle %d outside of grid (%g, %g), marking as inactive.\n",
                i, job->particles[i].x, job->particles[i].y);
            job->particles[i].active = 0;
        }


/*        if (job->in_element[i] < 0 || job->in_element[i] >= job->num_elements) {*/
/*            fprintf(stderr,*/
/*                "Particle %d outside of grid (%g, %g), marking as inactive.\n",*/
/*                i, job->particles[i].x, job->particles[i].y);*/
/*            job->particles[i].active = 0;*/
/*        }*/

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
        job->v_grid[2 * i] = job->nodes[i].x_t;
        job->v_grid[2 * i + 1] = job->nodes[i].y_t;
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void calculate_strainrate(job_t *job)
{
    int i;
    double dx_tdy;
    double dy_tdx;

    for (i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);

        job->particles[i].exx_t = DX_N_TO_P(job, x_t, 1.0, i);
        job->particles[i].eyy_t = DY_N_TO_P(job, y_t, 1.0, i);

        dx_tdy = DY_N_TO_P(job, x_t, 1.0, i);
        dy_tdx = DX_N_TO_P(job, y_t, 1.0, i);

        job->particles[i].exy_t = 0.5 * (dx_tdy + dy_tdx);
        job->particles[i].wxy_t = 0.5 * (dx_tdy - dy_tdx);
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void update_stress(job_t *job)
{
    int i;

    for (i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);

        ACCUMULATE_WITH_MUL(fx, job, sxx, i, n, b1, -job->particles[i].v);
        ACCUMULATE_WITH_MUL(fx, job, sxy, i, n, b2, -job->particles[i].v);
        ACCUMULATE_WITH_MUL(fy, job, sxy, i, n, b1, -job->particles[i].v);
        ACCUMULATE_WITH_MUL(fy, job, syy, i, n, b2, -job->particles[i].v);
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void update_internal_stress(job_t *job)
{
    int i;
    int p;
    int nn[4];

    double sxx, sxy, syy;
    double pxx, pxy, pyx, pyy;
    double dudx, dudy, dvdx, dvdy, jdet;

    for (i = 0; i < job->vec_len; i++) {
        job->f_int_grid[i] = 0;
    }

    for (i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);

        p = job->in_element[i];
        nn[0]  = job->elements[p].nodes[0];
        nn[1]  = job->elements[p].nodes[1];
        nn[2]  = job->elements[p].nodes[2];
        nn[3]  = job->elements[p].nodes[3];

        dudx = job->u_grid[2*nn[0]+0]*job->b11[i];
        dudx += job->u_grid[2*nn[1]+0]*job->b12[i];
        dudx += job->u_grid[2*nn[2]+0]*job->b13[i];
        dudx += job->u_grid[2*nn[3]+0]*job->b14[i];

        dudy = job->u_grid[2*nn[0]+0]*job->b21[i];
        dudy += job->u_grid[2*nn[1]+0]*job->b22[i];
        dudy += job->u_grid[2*nn[2]+0]*job->b23[i];
        dudy += job->u_grid[2*nn[3]+0]*job->b24[i];

        dvdx = job->u_grid[2*nn[0]+1]*job->b11[i];
        dvdx += job->u_grid[2*nn[1]+1]*job->b12[i];
        dvdx += job->u_grid[2*nn[2]+1]*job->b13[i];
        dvdx += job->u_grid[2*nn[3]+1]*job->b14[i];

        dvdy = job->u_grid[2*nn[0]+1]*job->b21[i];
        dvdy += job->u_grid[2*nn[1]+1]*job->b22[i];
        dvdy += job->u_grid[2*nn[2]+1]*job->b23[i];
        dvdy += job->u_grid[2*nn[3]+1]*job->b24[i];

        jdet = 1 + dudx + dvdy;

        sxx = job->particles[i].sxx;
        sxy = job->particles[i].sxy;
        syy = job->particles[i].syy;

        pxx = jdet * ((1 - dudx) * sxx - dudy * sxy);
        pxy = jdet * ((-dvdx) * sxx + (1 - dvdy) * sxy);
        pyx = jdet * ((1 - dudx) * sxy - dudy * syy);
        pyy = jdet * ((-dvdx) * sxy + (1 - dvdy) * syy);

        job->f_int_grid[2*nn[0]+0] += job->particles[i].v*(job->b11[i]*pxx + job->b21[i]*pxy);
        job->f_int_grid[2*nn[0]+1] += job->particles[i].v*(job->b11[i]*pyx + job->b21[i]*pyy);
        job->f_int_grid[2*nn[1]+0] += job->particles[i].v*(job->b12[i]*pxx + job->b22[i]*pxy);
        job->f_int_grid[2*nn[1]+1] += job->particles[i].v*(job->b12[i]*pyx + job->b22[i]*pyy);
        job->f_int_grid[2*nn[2]+0] += job->particles[i].v*(job->b13[i]*pxx + job->b23[i]*pxy);
        job->f_int_grid[2*nn[2]+1] += job->particles[i].v*(job->b13[i]*pyx + job->b23[i]*pyy);
        job->f_int_grid[2*nn[3]+0] += job->particles[i].v*(job->b14[i]*pxx + job->b24[i]*pxy);
        job->f_int_grid[2*nn[3]+1] += job->particles[i].v*(job->b14[i]*pyx + job->b24[i]*pyy);
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void build_stiffness_matrix(job_t *job)
{
    int i;
    int p;
    int nn[4];
    int lda = job->vec_len;
    double inv_dt_sq;

    /* zero matrix first */
    for (i = 0; i < (job->vec_len * job->vec_len); i++) {
        job->kku_grid[i] = 0;
    }

    /* build with K_geo and K_mat */
    for (i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);

        p = job->in_element[i];
        nn[0]  = job->elements[p].nodes[0];
        nn[1]  = job->elements[p].nodes[1];
        nn[2]  = job->elements[p].nodes[2];
        nn[3]  = job->elements[p].nodes[3];

    /* --- K_mat --- Symmetric? True --- */
    job->kku_grid[lda*(2*nn[0]+0)+2*nn[0]+0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*pow(job->b21[i], 2) - 1.0*pow(job->b11[i], 2) - 0.5*pow(job->b21[i], 2))/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[0]+0)+2*nn[0]+1] += ( 0.5*EMOD*job->b11[i]*job->b21[i]*job->particles[i].v/(-NUMOD + 1) );
    job->kku_grid[lda*(2*nn[0]+0)+2*nn[1]+0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b21[i]*job->b22[i] - 1.0*job->b11[i]*job->b12[i] - 0.5*job->b21[i]*job->b22[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[0]+0)+2*nn[1]+1] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b11[i]*job->b22[i] + 0.5*NUMOD*job->b12[i]*job->b21[i] - 0.5*job->b12[i]*job->b21[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[0]+0)+2*nn[2]+0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b21[i]*job->b23[i] - 1.0*job->b11[i]*job->b13[i] - 0.5*job->b21[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[0]+0)+2*nn[2]+1] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b11[i]*job->b23[i] + 0.5*NUMOD*job->b13[i]*job->b21[i] - 0.5*job->b13[i]*job->b21[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[0]+0)+2*nn[3]+0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b21[i]*job->b24[i] - 1.0*job->b11[i]*job->b14[i] - 0.5*job->b21[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[0]+0)+2*nn[3]+1] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b11[i]*job->b24[i] + 0.5*NUMOD*job->b14[i]*job->b21[i] - 0.5*job->b14[i]*job->b21[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[0]+1)+2*nn[0]+0] += ( 0.5*EMOD*job->b11[i]*job->b21[i]*job->particles[i].v/(-NUMOD + 1) );
    job->kku_grid[lda*(2*nn[0]+1)+2*nn[0]+1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*pow(job->b11[i], 2) - 0.5*pow(job->b11[i], 2) - 1.0*pow(job->b21[i], 2))/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[0]+1)+2*nn[1]+0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b22[i] - 1.0*NUMOD*job->b12[i]*job->b21[i] - 0.5*job->b11[i]*job->b22[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[0]+1)+2*nn[1]+1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b12[i] - 0.5*job->b11[i]*job->b12[i] - 1.0*job->b21[i]*job->b22[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[0]+1)+2*nn[2]+0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b23[i] - 1.0*NUMOD*job->b13[i]*job->b21[i] - 0.5*job->b11[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[0]+1)+2*nn[2]+1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b13[i] - 0.5*job->b11[i]*job->b13[i] - 1.0*job->b21[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[0]+1)+2*nn[3]+0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b24[i] - 1.0*NUMOD*job->b14[i]*job->b21[i] - 0.5*job->b11[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[0]+1)+2*nn[3]+1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b14[i] - 0.5*job->b11[i]*job->b14[i] - 1.0*job->b21[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[1]+0)+2*nn[0]+0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b21[i]*job->b22[i] - 1.0*job->b11[i]*job->b12[i] - 0.5*job->b21[i]*job->b22[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[1]+0)+2*nn[0]+1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b22[i] - 1.0*NUMOD*job->b12[i]*job->b21[i] - 0.5*job->b11[i]*job->b22[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[1]+0)+2*nn[1]+0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*pow(job->b22[i], 2) - 1.0*pow(job->b12[i], 2) - 0.5*pow(job->b22[i], 2))/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[1]+0)+2*nn[1]+1] += ( 0.5*EMOD*job->b12[i]*job->b22[i]*job->particles[i].v/(-NUMOD + 1) );
    job->kku_grid[lda*(2*nn[1]+0)+2*nn[2]+0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b22[i]*job->b23[i] - 1.0*job->b12[i]*job->b13[i] - 0.5*job->b22[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[1]+0)+2*nn[2]+1] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b12[i]*job->b23[i] + 0.5*NUMOD*job->b13[i]*job->b22[i] - 0.5*job->b13[i]*job->b22[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[1]+0)+2*nn[3]+0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b22[i]*job->b24[i] - 1.0*job->b12[i]*job->b14[i] - 0.5*job->b22[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[1]+0)+2*nn[3]+1] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b12[i]*job->b24[i] + 0.5*NUMOD*job->b14[i]*job->b22[i] - 0.5*job->b14[i]*job->b22[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[1]+1)+2*nn[0]+0] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b11[i]*job->b22[i] + 0.5*NUMOD*job->b12[i]*job->b21[i] - 0.5*job->b12[i]*job->b21[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[1]+1)+2*nn[0]+1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b12[i] - 0.5*job->b11[i]*job->b12[i] - 1.0*job->b21[i]*job->b22[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[1]+1)+2*nn[1]+0] += ( 0.5*EMOD*job->b12[i]*job->b22[i]*job->particles[i].v/(-NUMOD + 1) );
    job->kku_grid[lda*(2*nn[1]+1)+2*nn[1]+1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*pow(job->b12[i], 2) - 0.5*pow(job->b12[i], 2) - 1.0*pow(job->b22[i], 2))/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[1]+1)+2*nn[2]+0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b12[i]*job->b23[i] - 1.0*NUMOD*job->b13[i]*job->b22[i] - 0.5*job->b12[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[1]+1)+2*nn[2]+1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b12[i]*job->b13[i] - 0.5*job->b12[i]*job->b13[i] - 1.0*job->b22[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[1]+1)+2*nn[3]+0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b12[i]*job->b24[i] - 1.0*NUMOD*job->b14[i]*job->b22[i] - 0.5*job->b12[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[1]+1)+2*nn[3]+1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b12[i]*job->b14[i] - 0.5*job->b12[i]*job->b14[i] - 1.0*job->b22[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[2]+0)+2*nn[0]+0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b21[i]*job->b23[i] - 1.0*job->b11[i]*job->b13[i] - 0.5*job->b21[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[2]+0)+2*nn[0]+1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b23[i] - 1.0*NUMOD*job->b13[i]*job->b21[i] - 0.5*job->b11[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[2]+0)+2*nn[1]+0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b22[i]*job->b23[i] - 1.0*job->b12[i]*job->b13[i] - 0.5*job->b22[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[2]+0)+2*nn[1]+1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b12[i]*job->b23[i] - 1.0*NUMOD*job->b13[i]*job->b22[i] - 0.5*job->b12[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[2]+0)+2*nn[2]+0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*pow(job->b23[i], 2) - 1.0*pow(job->b13[i], 2) - 0.5*pow(job->b23[i], 2))/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[2]+0)+2*nn[2]+1] += ( 0.5*EMOD*job->b13[i]*job->b23[i]*job->particles[i].v/(-NUMOD + 1) );
    job->kku_grid[lda*(2*nn[2]+0)+2*nn[3]+0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b23[i]*job->b24[i] - 1.0*job->b13[i]*job->b14[i] - 0.5*job->b23[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[2]+0)+2*nn[3]+1] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b13[i]*job->b24[i] + 0.5*NUMOD*job->b14[i]*job->b23[i] - 0.5*job->b14[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[2]+1)+2*nn[0]+0] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b11[i]*job->b23[i] + 0.5*NUMOD*job->b13[i]*job->b21[i] - 0.5*job->b13[i]*job->b21[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[2]+1)+2*nn[0]+1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b13[i] - 0.5*job->b11[i]*job->b13[i] - 1.0*job->b21[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[2]+1)+2*nn[1]+0] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b12[i]*job->b23[i] + 0.5*NUMOD*job->b13[i]*job->b22[i] - 0.5*job->b13[i]*job->b22[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[2]+1)+2*nn[1]+1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b12[i]*job->b13[i] - 0.5*job->b12[i]*job->b13[i] - 1.0*job->b22[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[2]+1)+2*nn[2]+0] += ( 0.5*EMOD*job->b13[i]*job->b23[i]*job->particles[i].v/(-NUMOD + 1) );
    job->kku_grid[lda*(2*nn[2]+1)+2*nn[2]+1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*pow(job->b13[i], 2) - 0.5*pow(job->b13[i], 2) - 1.0*pow(job->b23[i], 2))/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[2]+1)+2*nn[3]+0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b13[i]*job->b24[i] - 1.0*NUMOD*job->b14[i]*job->b23[i] - 0.5*job->b13[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[2]+1)+2*nn[3]+1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b13[i]*job->b14[i] - 0.5*job->b13[i]*job->b14[i] - 1.0*job->b23[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[3]+0)+2*nn[0]+0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b21[i]*job->b24[i] - 1.0*job->b11[i]*job->b14[i] - 0.5*job->b21[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[3]+0)+2*nn[0]+1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b24[i] - 1.0*NUMOD*job->b14[i]*job->b21[i] - 0.5*job->b11[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[3]+0)+2*nn[1]+0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b22[i]*job->b24[i] - 1.0*job->b12[i]*job->b14[i] - 0.5*job->b22[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[3]+0)+2*nn[1]+1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b12[i]*job->b24[i] - 1.0*NUMOD*job->b14[i]*job->b22[i] - 0.5*job->b12[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[3]+0)+2*nn[2]+0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b23[i]*job->b24[i] - 1.0*job->b13[i]*job->b14[i] - 0.5*job->b23[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[3]+0)+2*nn[2]+1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b13[i]*job->b24[i] - 1.0*NUMOD*job->b14[i]*job->b23[i] - 0.5*job->b13[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[3]+0)+2*nn[3]+0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*pow(job->b24[i], 2) - 1.0*pow(job->b14[i], 2) - 0.5*pow(job->b24[i], 2))/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[3]+0)+2*nn[3]+1] += ( 0.5*EMOD*job->b14[i]*job->b24[i]*job->particles[i].v/(-NUMOD + 1) );
    job->kku_grid[lda*(2*nn[3]+1)+2*nn[0]+0] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b11[i]*job->b24[i] + 0.5*NUMOD*job->b14[i]*job->b21[i] - 0.5*job->b14[i]*job->b21[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[3]+1)+2*nn[0]+1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b14[i] - 0.5*job->b11[i]*job->b14[i] - 1.0*job->b21[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[3]+1)+2*nn[1]+0] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b12[i]*job->b24[i] + 0.5*NUMOD*job->b14[i]*job->b22[i] - 0.5*job->b14[i]*job->b22[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[3]+1)+2*nn[1]+1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b12[i]*job->b14[i] - 0.5*job->b12[i]*job->b14[i] - 1.0*job->b22[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[3]+1)+2*nn[2]+0] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b13[i]*job->b24[i] + 0.5*NUMOD*job->b14[i]*job->b23[i] - 0.5*job->b14[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[3]+1)+2*nn[2]+1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b13[i]*job->b14[i] - 0.5*job->b13[i]*job->b14[i] - 1.0*job->b23[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->kku_grid[lda*(2*nn[3]+1)+2*nn[3]+0] += ( 0.5*EMOD*job->b14[i]*job->b24[i]*job->particles[i].v/(-NUMOD + 1) );
    job->kku_grid[lda*(2*nn[3]+1)+2*nn[3]+1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*pow(job->b14[i], 2) - 0.5*pow(job->b14[i], 2) - 1.0*pow(job->b24[i], 2))/(pow(NUMOD, 2) - 1) );

    /* --- K_geo --- Symmetric? True --- */
    job->kku_grid[lda*(2*nn[0]+0)+2*nn[0]+0] += ( job->particles[i].v*(pow(job->b11[i], 2)*job->particles[i].sxx + 2*job->b11[i]*job->b21[i]*job->particles[i].sxy + pow(job->b21[i], 2)*job->particles[i].syy) );
    job->kku_grid[lda*(2*nn[0]+0)+2*nn[1]+0] += ( job->particles[i].v*(job->b11[i]*job->b12[i]*job->particles[i].sxx + job->b11[i]*job->b22[i]*job->particles[i].sxy + job->b12[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b22[i]*job->particles[i].syy) );
    job->kku_grid[lda*(2*nn[0]+0)+2*nn[2]+0] += ( job->particles[i].v*(job->b11[i]*job->b13[i]*job->particles[i].sxx + job->b11[i]*job->b23[i]*job->particles[i].sxy + job->b13[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b23[i]*job->particles[i].syy) );
    job->kku_grid[lda*(2*nn[0]+0)+2*nn[3]+0] += ( job->particles[i].v*(job->b11[i]*job->b14[i]*job->particles[i].sxx + job->b11[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b24[i]*job->particles[i].syy) );
    job->kku_grid[lda*(2*nn[0]+1)+2*nn[0]+1] += ( job->particles[i].v*(pow(job->b11[i], 2)*job->particles[i].sxx + 2*job->b11[i]*job->b21[i]*job->particles[i].sxy + pow(job->b21[i], 2)*job->particles[i].syy) );
    job->kku_grid[lda*(2*nn[0]+1)+2*nn[1]+1] += ( job->particles[i].v*(job->b11[i]*job->b12[i]*job->particles[i].sxx + job->b11[i]*job->b22[i]*job->particles[i].sxy + job->b12[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b22[i]*job->particles[i].syy) );
    job->kku_grid[lda*(2*nn[0]+1)+2*nn[2]+1] += ( job->particles[i].v*(job->b11[i]*job->b13[i]*job->particles[i].sxx + job->b11[i]*job->b23[i]*job->particles[i].sxy + job->b13[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b23[i]*job->particles[i].syy) );
    job->kku_grid[lda*(2*nn[0]+1)+2*nn[3]+1] += ( job->particles[i].v*(job->b11[i]*job->b14[i]*job->particles[i].sxx + job->b11[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b24[i]*job->particles[i].syy) );
    job->kku_grid[lda*(2*nn[1]+0)+2*nn[0]+0] += ( job->particles[i].v*(job->b11[i]*job->b12[i]*job->particles[i].sxx + job->b11[i]*job->b22[i]*job->particles[i].sxy + job->b12[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b22[i]*job->particles[i].syy) );
    job->kku_grid[lda*(2*nn[1]+0)+2*nn[1]+0] += ( job->particles[i].v*(pow(job->b12[i], 2)*job->particles[i].sxx + 2*job->b12[i]*job->b22[i]*job->particles[i].sxy + pow(job->b22[i], 2)*job->particles[i].syy) );
    job->kku_grid[lda*(2*nn[1]+0)+2*nn[2]+0] += ( job->particles[i].v*(job->b12[i]*job->b13[i]*job->particles[i].sxx + job->b12[i]*job->b23[i]*job->particles[i].sxy + job->b13[i]*job->b22[i]*job->particles[i].sxy + job->b22[i]*job->b23[i]*job->particles[i].syy) );
    job->kku_grid[lda*(2*nn[1]+0)+2*nn[3]+0] += ( job->particles[i].v*(job->b12[i]*job->b14[i]*job->particles[i].sxx + job->b12[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b22[i]*job->particles[i].sxy + job->b22[i]*job->b24[i]*job->particles[i].syy) );
    job->kku_grid[lda*(2*nn[1]+1)+2*nn[0]+1] += ( job->particles[i].v*(job->b11[i]*job->b12[i]*job->particles[i].sxx + job->b11[i]*job->b22[i]*job->particles[i].sxy + job->b12[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b22[i]*job->particles[i].syy) );
    job->kku_grid[lda*(2*nn[1]+1)+2*nn[1]+1] += ( job->particles[i].v*(pow(job->b12[i], 2)*job->particles[i].sxx + 2*job->b12[i]*job->b22[i]*job->particles[i].sxy + pow(job->b22[i], 2)*job->particles[i].syy) );
    job->kku_grid[lda*(2*nn[1]+1)+2*nn[2]+1] += ( job->particles[i].v*(job->b12[i]*job->b13[i]*job->particles[i].sxx + job->b12[i]*job->b23[i]*job->particles[i].sxy + job->b13[i]*job->b22[i]*job->particles[i].sxy + job->b22[i]*job->b23[i]*job->particles[i].syy) );
    job->kku_grid[lda*(2*nn[1]+1)+2*nn[3]+1] += ( job->particles[i].v*(job->b12[i]*job->b14[i]*job->particles[i].sxx + job->b12[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b22[i]*job->particles[i].sxy + job->b22[i]*job->b24[i]*job->particles[i].syy) );
    job->kku_grid[lda*(2*nn[2]+0)+2*nn[0]+0] += ( job->particles[i].v*(job->b11[i]*job->b13[i]*job->particles[i].sxx + job->b11[i]*job->b23[i]*job->particles[i].sxy + job->b13[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b23[i]*job->particles[i].syy) );
    job->kku_grid[lda*(2*nn[2]+0)+2*nn[1]+0] += ( job->particles[i].v*(job->b12[i]*job->b13[i]*job->particles[i].sxx + job->b12[i]*job->b23[i]*job->particles[i].sxy + job->b13[i]*job->b22[i]*job->particles[i].sxy + job->b22[i]*job->b23[i]*job->particles[i].syy) );
    job->kku_grid[lda*(2*nn[2]+0)+2*nn[2]+0] += ( job->particles[i].v*(pow(job->b13[i], 2)*job->particles[i].sxx + 2*job->b13[i]*job->b23[i]*job->particles[i].sxy + pow(job->b23[i], 2)*job->particles[i].syy) );
    job->kku_grid[lda*(2*nn[2]+0)+2*nn[3]+0] += ( job->particles[i].v*(job->b13[i]*job->b14[i]*job->particles[i].sxx + job->b13[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b23[i]*job->particles[i].sxy + job->b23[i]*job->b24[i]*job->particles[i].syy) );
    job->kku_grid[lda*(2*nn[2]+1)+2*nn[0]+1] += ( job->particles[i].v*(job->b11[i]*job->b13[i]*job->particles[i].sxx + job->b11[i]*job->b23[i]*job->particles[i].sxy + job->b13[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b23[i]*job->particles[i].syy) );
    job->kku_grid[lda*(2*nn[2]+1)+2*nn[1]+1] += ( job->particles[i].v*(job->b12[i]*job->b13[i]*job->particles[i].sxx + job->b12[i]*job->b23[i]*job->particles[i].sxy + job->b13[i]*job->b22[i]*job->particles[i].sxy + job->b22[i]*job->b23[i]*job->particles[i].syy) );
    job->kku_grid[lda*(2*nn[2]+1)+2*nn[2]+1] += ( job->particles[i].v*(pow(job->b13[i], 2)*job->particles[i].sxx + 2*job->b13[i]*job->b23[i]*job->particles[i].sxy + pow(job->b23[i], 2)*job->particles[i].syy) );
    job->kku_grid[lda*(2*nn[2]+1)+2*nn[3]+1] += ( job->particles[i].v*(job->b13[i]*job->b14[i]*job->particles[i].sxx + job->b13[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b23[i]*job->particles[i].sxy + job->b23[i]*job->b24[i]*job->particles[i].syy) );
    job->kku_grid[lda*(2*nn[3]+0)+2*nn[0]+0] += ( job->particles[i].v*(job->b11[i]*job->b14[i]*job->particles[i].sxx + job->b11[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b24[i]*job->particles[i].syy) );
    job->kku_grid[lda*(2*nn[3]+0)+2*nn[1]+0] += ( job->particles[i].v*(job->b12[i]*job->b14[i]*job->particles[i].sxx + job->b12[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b22[i]*job->particles[i].sxy + job->b22[i]*job->b24[i]*job->particles[i].syy) );
    job->kku_grid[lda*(2*nn[3]+0)+2*nn[2]+0] += ( job->particles[i].v*(job->b13[i]*job->b14[i]*job->particles[i].sxx + job->b13[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b23[i]*job->particles[i].sxy + job->b23[i]*job->b24[i]*job->particles[i].syy) );
    job->kku_grid[lda*(2*nn[3]+0)+2*nn[3]+0] += ( job->particles[i].v*(pow(job->b14[i], 2)*job->particles[i].sxx + 2*job->b14[i]*job->b24[i]*job->particles[i].sxy + pow(job->b24[i], 2)*job->particles[i].syy) );
    job->kku_grid[lda*(2*nn[3]+1)+2*nn[0]+1] += ( job->particles[i].v*(job->b11[i]*job->b14[i]*job->particles[i].sxx + job->b11[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b24[i]*job->particles[i].syy) );
    job->kku_grid[lda*(2*nn[3]+1)+2*nn[1]+1] += ( job->particles[i].v*(job->b12[i]*job->b14[i]*job->particles[i].sxx + job->b12[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b22[i]*job->particles[i].sxy + job->b22[i]*job->b24[i]*job->particles[i].syy) );
    job->kku_grid[lda*(2*nn[3]+1)+2*nn[2]+1] += ( job->particles[i].v*(job->b13[i]*job->b14[i]*job->particles[i].sxx + job->b13[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b23[i]*job->particles[i].sxy + job->b23[i]*job->b24[i]*job->particles[i].syy) );
    job->kku_grid[lda*(2*nn[3]+1)+2*nn[3]+1] += ( job->particles[i].v*(pow(job->b14[i], 2)*job->particles[i].sxx + 2*job->b14[i]*job->b24[i]*job->particles[i].sxy + pow(job->b24[i], 2)*job->particles[i].syy) );
    }

    inv_dt_sq = 1 / (job->dt * job->dt);

    /* diagonal mass matrix */
    for (i = 0; i < job->vec_len; i++) {
        if (job->m_grid[i] > TOL) {
            job->kku_grid[i + lda * i] += 4 * job->m_grid[i] * inv_dt_sq;
        } else {
            job->kku_grid[i + lda * i] += 1e-4;
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void implicit_solve(job_t *job)
{
    int i, j;

    /* lapack related*/
    int N;
    int lda, ldb;
    double *A;
    double *b;

    /* dnrm2 */
    int incx = 1;
    double du_norm = 0;
    double du0_norm = 0;
    double q_norm = 0;
    double q0_norm = 0;

    /* iteration count */
    int k = 0;
    static int stable_timestep = 0;

    /* old timestep grid velocity and acceleration */
    double *v_grid_t;
    double *a_grid_t;

    /* csparse */
    int nnz;
    int slda;
    cs *triplets;
    cs *smat;
    double *sb;

    /* Timing code */
    struct timespec requestStart, requestEnd;
    long ns;

    lda = job->vec_len;
    ldb = lda;
    N = lda;

    /* Lapack functions mutate matrices, so copy here first */
    A = (double *)malloc(lda * lda * sizeof(double));
    b = (double *)malloc(lda * sizeof(double));

    v_grid_t = (double *)malloc(lda * sizeof(double));
    a_grid_t = (double *)malloc(lda * sizeof(double));

    /* keep old timestep velocity and acceleration */
    memcpy(v_grid_t, job->v_grid, lda * sizeof(double));
    memcpy(a_grid_t, job->a_grid, lda * sizeof(double));

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

    /* Begin newton iterations */
    do {
        update_internal_stress(job);
        build_stiffness_matrix(job);

        /* reset stress state at beginning of newton iteration */
        for (i = 0; i < job->num_particles; i++) {
            CHECK_ACTIVE(job, i);
            
            memcpy(job->particles[i].state, job->particles[i].real_state, DEPVAR * sizeof(double));
            job->particles[i].sxx = job->particles[i].real_sxx;
            job->particles[i].sxy = job->particles[i].real_sxy;
            job->particles[i].syy = job->particles[i].real_syy;
        }

        /* Q = f_ext - f_int - M_g * (4 * u / dt^2 - 4 * v / dt - a) */
        memcpy(job->q_grid, job->m_grid, lda * sizeof(double));
        for (i = 0; i < lda; i++) {
            job->q_grid[i] *= -(4 * job->u_grid[i] / (job->dt * job->dt) - 4 * v_grid_t[i] / (job->dt) - a_grid_t[i]);
            job->q_grid[i] += (job->f_ext_grid[i] - job->f_int_grid[i]);
/*            job->q_grid[i] = (job->f_ext_grid[i] - job->f_int_grid[i]);*/
        }

#if 0
        for (i = 0; i < lda * lda; i++) {
            job->kku_grid[i] *= (job->dt * job->dt);
        }
        for (i = 0; i < lda; i++) {
            job->q_grid[i] *= (job->dt  * job->dt);
        }
#endif

        memcpy(A, job->kku_grid, lda * lda * sizeof(double));
        memcpy(b, job->q_grid, lda * sizeof(double));

        /* apply boundary conditions -- sticky + stiff bottom floor */
        for (i = 0; i < job->N; i++) {
            /* bottom floor */
            b[2 * i] = 0;
            b[2 * i + 1] = 0;
            for (j = 0; j < lda; j++) {
                A[lda * (j) + (2*i)] = 0;
                A[lda * (j) + (2*i + 1)] = 0;
            }
            A[lda * (2*i) + (2*i)] = 1;
            A[lda * (2*i + 1) + (2*i + 1)] = 1;

            /* left wall */
            b[2 * i * job->N] = 0;
/*            b[2 * i * job->N + 1] = 0;*/
            for (j = 0; j < lda; j++) {
                A[lda * (j) + (2*i*job->N)] = 0;
/*                A[lda * (j) + (2*i*job->N + 1)] = 0;*/
            }
            A[lda * (2*i*job->N) + (2*i*job->N)] = 1;
/*            A[lda * (2*i*job->N + 1) + (2*i*job->N + 1)] = 1;*/

            i++;
            /* right wall */
            b[2 * i * job->N - 2] = 0;
/*            b[2 * i * job->N - 1] = 0;*/
            for (j = 0; j < lda; j++) {
                A[lda * (j) + (2*i*job->N - 2)] = 0;
/*                A[lda * (j) + (2*i*job->N - 1)] = 0;*/
            }
            A[lda * (2*i*job->N - 2) + (2*i*job->N - 2)] = 1;
/*            A[lda * (2*i*job->N - 1) + (2*i*job->N - 1)] = 1;*/
            i--;

        }

        /* build node and inverse maps */
        for (i = 0; i < lda; i++) {
            job->node_u_map[i] = -1;
            job->inv_node_u_map[i] = -1;
        }
        j = 0;
        for (i = 0; i < lda; i++) {
            if (job->m_grid[i] > TOL) {
                job->node_u_map[i] = j;
                j++;
            }
        }
        slda = j;
        for (i = 0; i < lda; i++) {
            if ((j = job->node_u_map[i]) != -1) {
                job->inv_node_u_map[j] = i;
            }
        }

        /* build sparse matrix version */
        nnz = 0;
        for (i = 0; i < lda; i++) {
            for (j = 0; j < lda; j++) {
                if (abs(A[lda * i + j]) > TOL) {
                    nnz++;
                }
            }
        }

        triplets = cs_spalloc(slda, slda, nnz, 1, 1);
        sb = (double *)malloc(slda * sizeof(double));
        for (i = 0; i < slda; i++) {
            sb[i] = b[job->inv_node_u_map[i]];
        }

        int res = 0;
        for (j = 0; j < slda; j++) {
            for (i = 0; i < slda; i++) {
                if (abs(A[lda * job->inv_node_u_map[j] + job->inv_node_u_map[i]]) > TOL) {
                    res = cs_entry(triplets, i, j, A[lda * job->inv_node_u_map[j] + job->inv_node_u_map[i]]);
                    if (res == 0) {
                        fprintf(stderr, "error adding entry\n");
                    }
                }
            }
        }

/*        cs_print(triplets, 0);*/
        smat = cs_compress(triplets);
/*        cs_print(smat, 0);*/

        /* starting timer */
        clock_gettime(CLOCK_REALTIME, &requestStart);

        N = slda;
        if (k == 0) {
            q0_norm = dnrm2_(&N, sb, &incx);
        }
        q_norm = dnrm2_(&N, sb, &incx);

        if (!cs_lusol(0, smat, sb, 1e-12)) {
            fprintf(stderr, "lusol error!\n");
        }

        if (k == 0) {
            du0_norm = dnrm2_(&N, sb, &incx);
        }
        du_norm = dnrm2_(&N, sb, &incx);

        /* stopping timer */
        clock_gettime(CLOCK_REALTIME, &requestEnd);

        /* Calculate time it took */
    #define NS_PER_S 1E9
        ns = NS_PER_S * (requestEnd.tv_sec - requestStart.tv_sec )
          + ( requestEnd.tv_nsec - requestStart.tv_nsec );
        printf("Iteration[%d]: Norm du: %f Norm q: %f Implicit Solve: %lf s\n", k, du_norm, q_norm, (double)ns/NS_PER_S );
    #undef NS_PER_S

/*        memcpy(job->du_grid, b, lda * sizeof(double));*/
        for (i = 0; i < ldb; i++) {
            if ((j = job->node_u_map[i]) != -1) {
                job->du_grid[i] = sb[j];
            } else {
                job->du_grid[i] = 0;
            }
        }

        for (i = 0; i < ldb; i++) {
            /* update grid displacement */
            job->u_grid[i] += job->du_grid[i];

            /* update grid velocity */
            job->v_grid[i] = 2 * job->u_grid[i] / job->dt - v_grid_t[i];
        }

        for (i = 0; i < job->num_nodes; i++) {
            job->nodes[i].x_t = job->v_grid[2 * i];
            job->nodes[i].y_t = job->v_grid[2 * i + 1];
        }

        calculate_strainrate(job);
        calculate_stress(job);

        cs_spfree(triplets);
        cs_spfree(smat);
        free(sb);
        k++;
#define UNORM_TOL 1e-1
#define QNORM_TOL 1e-2
#define UNORM_CONV 1e-8

        if (du_norm < UNORM_CONV) {
            printf("norm of u  = %f: Accepting as converged.\n", du_norm);
            break;
        }

        if (q_norm*du_norm > 10*(q0_norm*du0_norm) || du_norm > 10*du0_norm || k > 10) {
            job->dt = job->dt * 0.8;
            stable_timestep = 0;
            printf("Trouble converging, modifying Timestep to %f.\n", job->dt);
            goto start_implicit;
        }

#define DT_LOWERBOUND 1e-8
        if (job->dt < DT_LOWERBOUND) {
            fprintf(stderr, "Timestep %g is too small, aborting.\n", job->dt);
            exit(-1);
        }

    } while ((du_norm / du0_norm) > UNORM_TOL || ((du_norm *q_norm) / (du0_norm * q0_norm)) > QNORM_TOL);

/*    stable_timestep++;*/
/*#define DT_UPPERBOUND 1.0/240.0*/
/*    if (stable_timestep > 4 && job->dt < DT_UPPERBOUND) {*/
/*        job->dt = job->dt * 1.25;*/
/*        if (job->dt > DT_UPPERBOUND) {*/
/*            job->dt = DT_UPPERBOUND;*/
/*        }*/
/*        printf("Increasing timestep to %f.\n", job->dt);*/
/*        stable_timestep = 0;*/
/*    }*/

    /* update grid acceleration */
    for (i = 0; i < ldb; i++) {
        job->a_grid[i] = (4 * job->u_grid[i] / (job->dt * job->dt) - 4 * v_grid_t[i] / (job->dt) - a_grid_t[i]);
    }

    /* repackage for use with macros */
    for (i = 0; i < job->num_nodes; i++) {
        job->nodes[i].x_tt = job->a_grid[2 * i];
        job->nodes[i].y_tt = job->a_grid[2 * i + 1];

        job->nodes[i].x_t = job->v_grid[2 * i];
        job->nodes[i].y_t = job->v_grid[2 * i + 1];

        job->nodes[i].ux = job->u_grid[2 * i];
        job->nodes[i].uy = job->u_grid[2 * i + 1];
    }

    free(A);
    free(b);
    free(v_grid_t);
    free(a_grid_t);

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
        job->m_grid[2 * i] = job->nodes[i].m;
        job->m_grid[2 * i + 1] = job->nodes[i].m;

        /* need previous timestep's acceleration for implicit method */
        if (job->nodes[i].m > TOL) {
            job->a_grid[2 * i] =  job->nodes[i].mx_tt / job->nodes[i].m;
            job->a_grid[2 * i + 1] = job->nodes[i].my_tt  / job->nodes[i].m;
        } else {
            job->a_grid[2 * i] =  0;
            job->a_grid[2 * i + 1] = 0;
        }

        /* external forces */
        job->f_ext_grid[2 * i] = job->nodes[i].fx;
        job->f_ext_grid[2 * i + 1] = job->nodes[i].fy;
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

        a_x_t = job->particles[i].x_tt;
        a_y_t = job->particles[i].y_tt;

        job->particles[i].x_tt = (N_TO_P(job, x_tt, i));
        job->particles[i].y_tt = (N_TO_P(job, y_tt, i));

        job->particles[i].x_t += 0.5 * job->dt * (job->particles[i].x_tt + a_x_t);
        job->particles[i].y_t += 0.5 * job->dt * (job->particles[i].y_tt + a_y_t);
/*        job->particles[i].x_t = 0.5 *(job->particles[i].x_t + (N_TO_P(job, x_t, i)));*/
/*        job->particles[i].y_t = 0.5 *(job->particles[i].y_t + (N_TO_P(job, y_t, i)));*/
    }
    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void update_deformation_gradient(job_t *job)
{
/*    #pragma omp parallel*/
/*    {*/
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

/*        #pragma omp for*/
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

            /* also update particle volume */
            job->particles[i].v = (job->particles[i].Fxx * job->particles[i].Fyy - job->particles[i].Fxy * job->particles[i].Fyx) * job->particles[i].v0;
        }
/*    }*/

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

        #define GRAD_THRESHOLD 0.5f

    for (i = 0; i < job->num_nodes; i++) {
        job->nodes[i].rho = job->nodes[i].mass_filled_element_neighbors / (job->nodes[i].num_filled_element_neighbors * job->h * job->h);
    }

    for (i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);
    /* XXX: Trololololol. Can't believe this looks so good :O
        -- not correct though, fix soon! */
/*        job->particles[i].v = job->h * job->h / (2 * job->elements[job->in_element[i]].n);*/

        /* Fixed version */
        job->particles[i].v = job->h * job->h / (job->elements[job->in_element[i]].n);
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
    free(job->h5);
    free(job->h6);
    free(job->h7);
    free(job->h8);
    free(job->h9);

    free(job->b11);
    free(job->b12);
    free(job->b13);
    free(job->b14);
    free(job->b15);
    free(job->b16);
    free(job->b17);
    free(job->b18);
    free(job->b19);

    free(job->b21);
    free(job->b22);
    free(job->b23);
    free(job->b24);
    free(job->b25);
    free(job->b26);
    free(job->b27);
    free(job->b28);
    free(job->b29);

    free(job);

    return;
}
/*----------------------------------------------------------------------------*/

