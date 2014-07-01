/**
    \file g_local.c
    \author Sachith Dunatunga
    \date 23.10.13

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

#define signum(x) ((int)((0 < x) - (x < 0)))

#define jp(x) job->particles[i].x

#define mu0     0.2f
#define hmu     360.0f
#define p       1.88f
#define mucv    0.613f
#define etacr   0.54f
#define b       2.25f
#define q       1.0f
#define eta0    0.675f
#define hbeta   1.5f

#undef EMOD
#undef NUMOD

#define EMOD 1e8
#define NUMOD 0.3

#define G (EMOD / (2.0f * (1.0f + NUMOD)))
#define K (EMOD / (3.0f * (1.0f - 2*NUMOD)))

#define PI 3.1415926535897932384626433
#define PHI (30.0*PI/180.0)

/* from Dave + Ken's paper (modified B) */
#define MU_S 0.3819
#define GRAINS_RHO 2450
#define B (0.9377 * 20)

/*
    from geometric considerations -- artificially increased to make
    velocity field reasonable in chute.
*/
#define GRAINS_D (0.01 * 5)

#define Epxx jp(state[0])
#define Epxy jp(state[1])
#define Epyy jp(state[2])
#define gf jp(state[3])
#define eta jp(state[4])
#define beta jp(state[5])
#define gammap jp(state[9])
#define Ef_mag jp(state[10])
#define Exx jp(state[6])
#define Exy jp(state[7])
#define Eyy jp(state[8])

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
        Epxx = 0;
        Epxy = 0;
        Epyy = 0;
        Exx = 0;
        Exy = 0;
        Eyy = 0;
        eta = 0;
        gammap = 0;
        gf = 0;
    }

    printf("%s:%s: done initializing material.\n", __FILE__,  __func__);
    return;
}
/*----------------------------------------------------------------------------*/

/* Local granular fluidity model. */
/*----------------------------------------------------------------------------*/
void calculate_stress(job_t *job)
{
    double t0xx;
    double t0xy;
    double t0yy;

    double p_t;
    double tau_t;

    double dsjxx;
    double dsjxy;
    double dsjyy;


    double Np_trxx;
    double Np_trxy;
    double Np_tryy;

    double nup_tau;

    double dpxx;
    double dpxy;
    double dpyy;

    double const c = 1e-4;
    int i;

    for (i = 0; i < job->num_particles; i++) {
        if (job->active[i] == 0) {
            continue;
        }

        /* check if the density allows for supporting any stress */
        if ((job->particles[i].m / job->particles[i].v) < 1485.0f) {
            job->particles[i].sxx = 0;
            job->particles[i].sxy = 0;
            job->particles[i].syy = 0;
            continue;
        }

        /* Calculate tau, p, and mu at beginning of timestep. */
        p_t = -0.5 * (job->particles[i].sxx + job->particles[i].syy);
        t0xx = job->particles[i].sxx + p_t;
        t0xy = job->particles[i].sxy;
        t0yy = job->particles[i].syy + p_t;
        tau_t = sqrt(0.5f*(t0xx*t0xx + 2*t0xy*t0xy + t0yy*t0yy));
        if (p_t > c) {
            eta = B * GRAINS_D * sqrt (GRAINS_RHO * p_t);
        } else {
            eta = B * GRAINS_D * sqrt (GRAINS_RHO * c);
        }

        /* Calculate direction of plastic flow. */
        if (tau_t > 0) {
            Np_trxx = sqrt(0.5) * t0xx / tau_t;
            Np_trxy = sqrt(0.5) * t0xy / tau_t;
            Np_tryy = sqrt(0.5) * t0yy / tau_t;
        } else {
            Np_trxx = 0;
            Np_trxy = 0;
            Np_tryy = 0;
        }
 
        if (tau_t > (MU_S * p_t) && p_t > c) {
            nup_tau = ((tau_t - MU_S * p_t) / eta);
            dpxx = sqrt(0.5) * nup_tau * Np_trxx;
            dpxy = sqrt(0.5) * nup_tau * Np_trxy;
            dpyy = sqrt(0.5) * nup_tau * Np_tryy;
        } else {
            dpxx = 0;
            dpxy = 0;
            dpyy = 0;
            nup_tau = 0;
        }

/*        if (p_t < c) {*/
/*            dpxx = job->particles[i].exx_t;*/
/*            dpxy = job->particles[i].exy_t;*/
/*            dpyy = job->particles[i].eyy_t;*/
/*            job->particles[i].sxx = -c;*/
/*            job->particles[i].sxy = 0;*/
/*            job->particles[i].syy = -c;*/
/*            continue;*/
/*        }*/

        /* use strain rate to calculate stress increment */
        dsjxx = job->dt * (EMOD / (1 - NUMOD*NUMOD)) * ((job->particles[i].exx_t - dpxx) + NUMOD * (job->particles[i].eyy_t - dpyy));
        dsjxy = job->dt * (EMOD / (2 *(1 + NUMOD))) * (job->particles[i].exy_t - dpxy);
        dsjyy = job->dt * (EMOD / (1 - NUMOD*NUMOD)) * ((job->particles[i].eyy_t - dpyy) + NUMOD * (job->particles[i].exx_t - dpxx));
        dsjxx -= 2 * job->dt * job->particles[i].wxy_t * job->particles[i].sxy;
        dsjxy += job->dt * job->particles[i].wxy_t * (job->particles[i].sxx - job->particles[i].syy);
        dsjyy += 2 * job->dt * job->particles[i].wxy_t * job->particles[i].sxy;

        job->particles[i].sxx += dsjxx;
        job->particles[i].sxy += dsjxy;
        job->particles[i].syy += dsjyy;

        gammap += nup_tau * job->dt;
    }

    return;
}
/*----------------------------------------------------------------------------*/

