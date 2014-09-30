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

#include "tensor.h"

#define jp(x) job->particles[i].x

#undef EMOD
#undef NUMOD

#undef G
#undef K

/* from Dave + Ken's paper */
#define MU_S 0.3819
#define GRAINS_RHO 2450

/* from Jop (modified I_0) */
#define MU_2 0.6435
#define I_0 (0.278 / 1.0)

/*
    from geometric considerations -- artificially increased to make
    velocity field reasonable in chute.
*/
#define GRAINS_D (0.001 * 5)

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
            /* what */
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
        x = (-b - sqrt(b*b - 4*a*c)) / (2*a);
    } else {
        x = (2*c) / (-b + sqrt(b*b - 4*a*c));
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
        job->particles[i].T[XX] = job->particles[i].sxx;
        job->particles[i].T[XY] = job->particles[i].sxy;
        job->particles[i].T[XZ] = 0;
        job->particles[i].T[YX] = job->particles[i].sxy;
        job->particles[i].T[YY] = job->particles[i].syy;
        job->particles[i].T[YZ] = 0;
        job->particles[i].T[ZX] = 0;
        job->particles[i].T[ZY] = 0;
        job->particles[i].T[ZZ] = 0;
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

    int density_flag;
    
    size_t i;

    double trD;
    const double lambda = K - 2.0 * G / 3.0;

    double S0, S2;
    double B, H;
    double alpha;

/*    fprintf(stderr, "processing particle ids [%zu %zu].\n", p_start, p_stop);*/

    for (i = p_start; i < p_stop; i++) {
        if (job->active[i] == 0) {
            continue;
        }

        double T0[9], Ttr[9]; /* deviator and trial */
        double JS[9]; /* Jaumann spin term. */
        double W[9], D[9]; /* spin and stretching */
        double temp[9];

        /* 3D velocity gradient (plane strain). */
        job->particles[i].L[XX] = job->particles[i].exx_t;
        job->particles[i].L[XY] = job->particles[i].exy_t + job->particles[i].wxy_t;
        job->particles[i].L[YX] = job->particles[i].exy_t - job->particles[i].wxy_t;
        job->particles[i].L[YY] = job->particles[i].eyy_t;
        job->particles[i].L[XZ] = 0;
        job->particles[i].L[YZ] = 0;
        job->particles[i].L[ZX] = 0;
        job->particles[i].L[ZY] = 0;
        job->particles[i].L[ZZ] = 0;

        tensor_skw(W, job->particles[i].L, 3);
        tensor_sym(D, job->particles[i].L, 3);
        tensor_trace(&trD, D, 3);

        /* construct jaumman spin term */        
        tensor_multiply(JS, W, job->particles[i].T, 3);
        tensor_multiply(temp, job->particles[i].T, W, 3);
        tensor_scale(temp, -1.0, 3);
        tensor_add(JS, JS, temp, 3);

        /* construct trial stress */
        tensor_copy(Ttr, D, 3);
        tensor_scale(Ttr, 2.0 * G, 3);
        Ttr[XX] += lambda * trD;
        Ttr[YY] += lambda * trD;
        Ttr[ZZ] += lambda * trD;
        tensor_add(Ttr, Ttr, JS, 3);
        tensor_scale(Ttr, job->dt, 3);
        tensor_add(Ttr, Ttr, job->particles[i].T, 3);

        /* Calculate tau and p trial values. */
        tensor_decompose(T0, &p_tr, Ttr, 3);
        p_tr *= -1.0;
        tensor_contraction(&tau_tr, T0, T0, 3);
        tau_tr *= 0.5;
        tau_tr = sqrt(tau_tr);

#if 0
        trD = job->particles[i].exx_t + job->particles[i].eyy_t;
        dsjxx = lambda * trD + 2.0 * G * job->particles[i].exx_t;
        dsjxy = 2.0 * G * job->particles[i].exy_t;
        dsjyy = lambda * trD + 2.0 * G * job->particles[i].eyy_t;
        dsjxx += 2 * job->particles[i].wxy_t * job->particles[i].sxy;
        dsjxy -= job->particles[i].wxy_t * (job->particles[i].sxx - job->particles[i].syy);
        dsjyy -= 2 * job->particles[i].wxy_t * job->particles[i].sxy;

        sxx_tr = job->particles[i].sxx + job->dt * dsjxx;
        sxy_tr = job->particles[i].sxy + job->dt * dsjxy;
        syy_tr = job->particles[i].syy + job->dt * dsjyy;

        p_tr = -0.5 * (sxx_tr + syy_tr);
        t0xx_tr = sxx_tr + p_tr;
        t0xy_tr = sxy_tr;
        t0yy_tr = syy_tr + p_tr;
        tau_tr = sqrt(0.5*(t0xx_tr*t0xx_tr + 2*t0xy_tr*t0xy_tr + t0yy_tr*t0yy_tr));
#endif

        if ((job->particles[i].m / job->particles[i].v) < 1485.0f) {
            density_flag = 1;
/*            printf("%4d: density %lf\n", i, (job->particles[i].m / job->particles[i].v));*/
        } else {
            density_flag = 0;
        }

        if (density_flag || p_tr <= c) {
            nup_tau = (tau_tr) / (G * job->dt);
            
            /* job->particles[i].sxx = 0;
            job->particles[i].sxy = 0;
            job->particles[i].syy = 0; */
            job->particles[i].T[XX] = 0;
            job->particles[i].T[XY] = 0;
            job->particles[i].T[XZ] = 0;
            job->particles[i].T[YX] = 0;
            job->particles[i].T[YY] = 0;
            job->particles[i].T[YZ] = 0;
            job->particles[i].T[ZX] = 0;
            job->particles[i].T[ZY] = 0;
            job->particles[i].T[ZZ] = 0;
        } else if (p_tr > c) {
            S0 = MU_S * p_tr;
            if (tau_tr <= S0) {
                tau_tau = tau_tr;
                scale_factor = 1.0;
            } else {
                S2 = MU_2 * p_tr;
                alpha = G * I_0 * job->dt * sqrt(p_tr / GRAINS_RHO) / GRAINS_D;
                B = -(S2 + tau_tr + alpha);
                H = S2 * tau_tr + S0 * alpha;
                tau_tau = negative_root(1.0, B, H);
                scale_factor = (tau_tau / tau_tr);
            }

            nup_tau = ((tau_tr - tau_tau) / G) / job->dt;

            tensor_scale(T0, scale_factor, 3);
            tensor_copy(job->particles[i].T, T0, 3);
            job->particles[i].T[XX] -= p_tr;
            job->particles[i].T[YY] -= p_tr;
            job->particles[i].T[ZZ] -= p_tr;
            /* job->particles[i].sxx = scale_factor * t0xx_tr - p_tr;
            job->particles[i].sxy = scale_factor * t0xy_tr;
            job->particles[i].syy = scale_factor * t0yy_tr - p_tr; */
        } else {
/*            fprintf(stderr, "u %zu %3.3g %3.3g %d ", i, f, p_tr, density_flag);*/
            fprintf(stderr, "u"); 
            nup_tau = 0;
        }

        /* Copy relevant stress entries. */
        job->particles[i].sxx = job->particles[i].T[XX];
        job->particles[i].sxy = job->particles[i].T[XY];
        job->particles[i].syy = job->particles[i].T[YY];
        
        /* use strain rate to calculate stress increment */
        gammap += nup_tau * job->dt;
        gammadotp = nup_tau;
    }

    return;
}
/*----------------------------------------------------------------------------*/

