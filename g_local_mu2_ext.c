/**
    \file g_local_mu2.c
    \author Sachith Dunatunga
    \date 04.12.13

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
#include "exitcodes.h"

#include <assert.h>

#define signum(x) ((int)((0 < x) - (x < 0)))

#define jp(x) job->particles[i].x

#define PI 3.1415926535897932384626433
#define PHI (30.0*PI/180.0)

/* from Dave + Ken's paper (modified B) */
#define MU_S 0.3819
#define GRAINS_RHO 2450
/*#define B (0.9377 * 20)*/
#define A 4.8

/* from Jop (modified I_0) */
#define MU_2 0.6435
#define I_0 (0.278 / 1.0)

/*#define MU_2 (30.0*PI/180.0)*/
/*#define MU_S MU_2*/

/*
    from geometric considerations -- artificially increased to make
    velocity field reasonable in chute.
*/
#define GRAINS_D (0.01 * 5)

#define mu_t jp(state[0])
/*#define Epxy jp(state[1])*/
/*#define Epyy jp(state[2])*/
#define gf jp(state[3])
#define gflocal jp(state[4])
#define beta jp(state[5])
#define gammap jp(state[9])
#define gammadotp jp(state[10])
#define sxx_e jp(state[6])
#define sxy_e jp(state[7])
#define syy_e jp(state[8])

#define GFLOCAL_IDX 4

#define MAT_VERSION_STRING "1.0 " __DATE__ " " __TIME__

void calculate_stress(job_t *job);
void calculate_stress_threaded(threadtask_t *task);

/*
    The Young's Modulus (E) and Poisson ratio (nu), set in the material init
    procedure. These are Ccpies of the properties given in the
    configuration file. Shear modulus (G) and bulk modulus (K) are derived from
    these.
*/
static double E;
static double nu;
static double G;
static double K;

/* granular fluidity at nodes. */
static double *gf_nodes;

/* laplacian of granular fluidity at nodes. */
static double *d2_gf_nodes;

inline double mu_from_gammadot(double gammadot, double p,
    double dsqrtrhos, double mu_s, double mu_2, double inum_0);

inline double mu_from_gammadot(double gammadot, double p,
    double dsqrtrhos, double mu_s, double mu_2, double inum_0)
{
    double mu;
    double inum;

    if (p <= 0) {
        mu = mu_2;
    } else {
        inum = dsqrtrhos * gammadot / sqrt(p);
        if (inum <= 0) {
            mu = mu_s;
        } else {
            mu = mu_s + (mu_2 - mu_s) / ((inum_0 / inum) + 1.0);
        }
    }

    return mu;
}

/*----------------------------------------------------------------------------*/
void material_init(job_t *job)
{
    size_t i, j;

    for (i = 0; i < job->num_particles; i++) {
        for (j = 0; j < DEPVAR; j++) {
            job->particles[i].state[j] = 0;
        }
    }

    for (i = 0; i < job->num_particles; i++) {
        mu_t = 0;
        sxx_e = 0;
        sxy_e = 0;
        syy_e = 0;
        beta = 0;
        gammap = 0;
        gammadotp = 0;
        gf = 0;
        gflocal = 0;
    }

    if (job->material.num_fp64_props < 2) {
        fprintf(stderr,
            "%s:%s: Need at least 2 properties defined (E, nu).\n",
            __FILE__, __func__);
        exit(EXIT_ERROR_MATERIAL_FILE);
    } else {
        E = job->material.fp64_props[0];
        nu = job->material.fp64_props[1];
        G = E / (2.0 * (1.0 + nu));
        K = E / (3.0 * (1.0 - 2*nu));
        printf("%s:%s: properties (E = %g, nu = %g, G = %g, K = %g).\n",
            __FILE__, __func__, E, nu, G , K);
    }

    gf_nodes = (double *)malloc(job->num_nodes * sizeof(double));
    d2_gf_nodes = (double *)malloc(job->num_nodes * sizeof(double));
    assert(gf_nodes != NULL);
    assert(d2_gf_nodes != NULL);

    printf("%s:%s: (material version %s) done initializing material.\n",
        __FILE__,  __func__, MAT_VERSION_STRING);
    return;
}
/*----------------------------------------------------------------------------*/

/* Local granular fluidity model. */
void calculate_stress(job_t *job)
{
    threadtask_t t;
    t.job = job;
    t.offset = 0;
    t.blocksize = job->num_particles;
    calculate_stress_threaded(&t);
    return;
}

/*----------------------------------------------------------------------------*/
void calculate_stress_threaded(threadtask_t *task)
{
    job_t *job = task->job;

    /* Since this is local, we can split the particles among the threads. */
    size_t p_start = task->offset;
    size_t p_stop = task->offset + task->blocksize;

    size_t n_start = task->n_offset;
    size_t n_stop = task->n_offset + task->n_blocksize;

    /* values from previous timestep */
    double p_t;
    /* kept as a state variable */
/*    double mu_t; */

    /* value at end of timestep */
    double tau_tau;
    double scale_factor;

    /* increment using jaumann rate */
    double dsjxx, dsjxy, dsjyy;

    /* trial values */
    double sxx_tr, sxy_tr, syy_tr;
    double t0xx_tr, t0xy_tr, t0yy_tr;
    double p_tr, tau_tr;

    double nup_tau;

    double const c = 1e-2;

    double f;
    int density_flag;
    
    size_t i, j, p;
    size_t cc, p_idx, tc_idx;

    double s[4];

/*    fprintf(stderr, "processing particle ids [%zu %zu].\n", p_start, p_stop);*/

    /* clear node-level values */
    for (i = n_start; i < n_stop; i++) {
        d2_gf_nodes[i] = 0;
    }

    for (i = p_start; i < p_stop; i++) {
        if (job->active[i] == 0) {
            continue;
        }

        /* Calculate p at beginning of timestep. */
        p_t = -0.5 * (job->particles[i].sxx + job->particles[i].syy);

        if (gammadotp == 0) {
            gflocal = 0;
        } else {
            mu_t = mu_from_gammadot(gammadotp, p_t,
                GRAINS_D * sqrt(GRAINS_RHO), MU_S, MU_2, I_0);
            
            gflocal = gammadotp / mu_t;
        }
    }
    pthread_barrier_wait(job->serialize_barrier);

    for (cc = 0; cc < job->num_colors; cc++) {
        tc_idx = task->id * job->num_colors + cc;
        for (i = 0; i < job->particle_by_element_color_lengths[tc_idx]; i++) {
            p_idx = job->particle_by_element_color_lists[tc_idx][i];

            /*
                Note: We don't need to check if the particle is active since
                the list assembly already checks for this condition.
            */

            p = job->in_element[p_idx];

            s[0] = job->b11[p_idx] * job->b11[p_idx] + job->b21[p_idx] * job->b21[p_idx];
            s[1] = job->b12[p_idx] * job->b12[p_idx] + job->b22[p_idx] * job->b22[p_idx];
            s[2] = job->b13[p_idx] * job->b13[p_idx] + job->b23[p_idx] * job->b23[p_idx];
            s[3] = job->b14[p_idx] * job->b14[p_idx] + job->b24[p_idx] * job->b24[p_idx];

            for (j = 0; j < 4; j++) {
                d2_gf_nodes[job->elements[p].nodes[j]] +=
                    -job->particles[p_idx].v * job->particles[p_idx].state[GFLOCAL_IDX] * s[j];
            }
        }

        /*
            Each color can be done simultaneously, but we have to sync between
            colors.
        */
        pthread_barrier_wait(job->serialize_barrier);
    }
/*    pthread_barrier_wait(job->serialize_barrier);*/

    for (i = p_start; i < p_stop; i++) {
        if (job->active[i] == 0) {
            continue;
        }

        gf = 0;

        p = job->in_element[i];
        for (j = 0; j < 4; j++) {
            gf += 0.25 * d2_gf_nodes[job->elements[p].nodes[j]];
        }

        gf = A * A * GRAINS_D * GRAINS_D * mu_t * gf;

        if (mu_t > MU_S) {
            gammadotp = gammadotp - (gf / fabs(mu_t - MU_S));
        }

        if (gammadotp < 0) {
            gammadotp = 0;
        }
    }


    pthread_barrier_wait(job->serialize_barrier);

    for (i = p_start; i < p_stop; i++) {
        if (job->active[i] == 0) {
            continue;
        }

        /* Calculate p at beginning of timestep. */
        p_t = -0.5 * (job->particles[i].sxx + job->particles[i].syy);

        /* Calculate tau and p trial values. */
        dsjxx = job->dt * (E / (1 - nu*nu)) * ((job->particles[i].exx_t) + nu * (job->particles[i].eyy_t));
        dsjxy = job->dt * (E / (2 *(1 + nu))) * (job->particles[i].exy_t);
        dsjyy = job->dt * (E / (1 - nu*nu)) * ((job->particles[i].eyy_t) + nu * (job->particles[i].exx_t));
        dsjxx -= 2 * job->dt * job->particles[i].wxy_t * job->particles[i].sxy;
        dsjxy += job->dt * job->particles[i].wxy_t * (job->particles[i].sxx - job->particles[i].syy);
        dsjyy += 2 * job->dt * job->particles[i].wxy_t * job->particles[i].sxy;

        sxx_tr = job->particles[i].sxx + dsjxx;
        sxy_tr = job->particles[i].sxy + dsjxy;
        syy_tr = job->particles[i].syy + dsjyy;

        p_tr = -0.5 * (sxx_tr + syy_tr);
        t0xx_tr = sxx_tr + p_tr;
        t0xy_tr = sxy_tr;
        t0yy_tr = syy_tr + p_tr;
        tau_tr = sqrt(0.5*(t0xx_tr*t0xx_tr + 2*t0xy_tr*t0xy_tr + t0yy_tr*t0yy_tr));

        mu_t = mu_from_gammadot(gammadotp, p_t,
                GRAINS_D * sqrt(GRAINS_RHO), MU_S, MU_2, I_0);

        tau_tau = mu_t*(p_tr + c);
        f = tau_tr - tau_tau;

        if ((job->particles[i].m / job->particles[i].v) < 1485.0f) {
            density_flag = 1;
/*            printf("%4d: density %lf\n", i, (job->particles[i].m / job->particles[i].v));*/
        } else {
            density_flag = 0;
        }

        if (density_flag) {
            nup_tau = (tau_tr) / (G * job->dt);
            beta = -p_tr / (K * job->dt);

            job->particles[i].sxx = 0;
            job->particles[i].sxy = 0;
            job->particles[i].syy = 0;
        } else if (f < 0) {
            nup_tau = 0;
            beta = 0;

            job->particles[i].sxx = t0xx_tr - p_tr;
            job->particles[i].sxy = t0xy_tr;
            job->particles[i].syy = t0yy_tr - p_tr;
        } else if (f >= 0 && p_tr > -c/mu_t) {
            nup_tau = (tau_tr - tau_tau) / (G * job->dt);
            beta = 0;

            scale_factor = (tau_tau / tau_tr);
            job->particles[i].sxx = scale_factor * t0xx_tr - p_tr;
            job->particles[i].sxy = scale_factor * t0xy_tr;
            job->particles[i].syy = scale_factor * t0yy_tr - p_tr;
        } else if (p_tr <= -c/mu_t) {
            nup_tau = tau_tr / (G * job->dt);
            beta = ((c/mu_t) - p_tr) / (K * job->dt);

            job->particles[i].sxx = 0;
            job->particles[i].sxy = 0;
            job->particles[i].syy = 0;
        } else {
            fprintf(stderr, "u %zu %3.3g %3.3g %d ", i, f, p_tr, density_flag);
/*            fprintf(stderr, "u"); */
            nup_tau = 0;
        }

        /* use strain rate to calculate stress increment */
        gammap += nup_tau * job->dt;
        gammadotp = nup_tau;
    }

    return;
}
/*----------------------------------------------------------------------------*/

