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

#define signum(x) ((int)((0 < x) - (x < 0)))

#define jp(x) job->particles[i].x

#undef EMOD
#undef NUMOD

#define G (E / (2.0f * (1.0f + nu)))
#define K (E / (3.0f * (1.0f - 2*nu)))

#undef G
#undef K

#define PI 3.1415926535897932384626433
#define PHI (30.0*PI/180.0)

/* from Dave + Ken's paper (modified B) */
#define MU_S 0.3819
#define GRAINS_RHO 2450
#define B (0.9377 * 20)

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

#define mu_y jp(state[0])
/*#define Epxy jp(state[1])*/
/*#define Epyy jp(state[2])*/
#define gf jp(state[3])
#define eta jp(state[4])
#define beta jp(state[5])
#define gammap jp(state[9])
#define gammadotp jp(state[10])
#define sxx_e jp(state[6])
#define sxy_e jp(state[7])
#define syy_e jp(state[8])

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

void quadratic_roots(double *x1, double *x2, double a, double b, double c)
{
    if (a == 0) {
        /* not a quadratic... */
        if (b == 0) {
            /* wut */
            *x1 = 0;
            *x2 = 0;
        } else {
            *x1 = c / b;
            *x2 = 0;
        }
    } else {
        *x1 = (-b - copysign(sqrt(b * b - 4 * a * c), b)) / (2 * a);
        *x2 = c / (a * *x1);
    }
    return;
}

double negative_root(double a, double b, double c)
{
    double x;
    if (b > 0) {
        x = (-b - sqrt(b*b - 4*a*c)) / (2 * a);
    } else {
        x = (2 * c) / (-b + sqrt(b*b - 4*a*c));
    }
    return x;
}

/*----------------------------------------------------------------------------*/
void material_init(job_t *job)
{
    int i, j;

    for (i = 0; i < job->num_particles; i++) {
        for (j = 0; j < DEPVAR; j++) {
            job->particles[i].state[j] = 0;
        }
    }

    for (i = 0; i < job->num_particles; i++) {
        mu_y = 0;
        sxx_e = 0;
        sxy_e = 0;
        syy_e = 0;
        eta = 0;
        beta = 0;
        gammap = 0;
        gammadotp = 0;
        gf = 0;
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

    /* values from previous timestep */
    double p_t, mu_t;

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

    double const c = 0;

    double f;
    int density_flag;
    
    size_t i;

    double inertial_num;
    double trD;
    const double lambda = K - 2.0 * G / 3.0;


/*    fprintf(stderr, "processing particle ids [%zu %zu].\n", p_start, p_stop);*/

    for (i = p_start; i < p_stop; i++) {
        if (job->active[i] == 0) {
            continue;
        }

        /* Calculate p at beginning of timestep. */
        p_t = -0.5 * (job->particles[i].sxx + job->particles[i].syy);

        /* Calculate tau and p trial values. */
        trD = job->particles[i].exx_t + job->particles[i].eyy_t;
        dsjxx = lambda * trD + 2.0 * G * job->particles[i].exx_t;
        dsjxy = 2.0 * G * job->particles[i].exy_t;
        dsjyy = lambda * trD + 2.0 * G * job->particles[i].eyy_t;
        dsjxx -= 2 * job->particles[i].wxy_t * job->particles[i].sxy;
        dsjxy += job->particles[i].wxy_t * (job->particles[i].sxx - job->particles[i].syy);
        dsjyy += 2 * job->particles[i].wxy_t * job->particles[i].sxy;

        sxx_tr = job->particles[i].sxx + job->dt * dsjxx;
        sxy_tr = job->particles[i].sxy + job->dt * dsjxy;
        syy_tr = job->particles[i].syy + job->dt * dsjyy;

        p_tr = -0.5 * (sxx_tr + syy_tr);
        t0xx_tr = sxx_tr + p_tr;
        t0xy_tr = sxy_tr;
        t0yy_tr = syy_tr + p_tr;
        tau_tr = sqrt(0.5*(t0xx_tr*t0xx_tr + 2*t0xy_tr*t0xy_tr + t0yy_tr*t0yy_tr));

        inertial_num = gammadotp * GRAINS_D * sqrt(GRAINS_RHO);
        if (p_t <= 0) {
            mu_t = MU_2;
        } else {
            inertial_num = inertial_num / sqrt(p_t);
            if (inertial_num <= 0) {
                mu_t = MU_S;
            } else {
                mu_t = MU_S + (MU_2 - MU_S) / ((I_0 / inertial_num) + 1.0);
            }
        }

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
/*            nup_tau = (tau_tr - tau_tau) / (G * job->dt);*/
/*            nup_tau = (tau_tr / (G * job->dt)) - (tau_tau / (G * job->dt));*/
            nup_tau = (f / G) / job->dt;
            beta = 0;

/*            printf("%g %g\n", f, nup_tau);*/

/*            if (nup_tau < 1e-8) {*/
/*                printf("p %zu %g %g %g ", i, nup_tau, tau_tau, tau_tr);*/
/*            }*/

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

