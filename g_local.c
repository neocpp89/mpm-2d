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
    double txx;
    double txy;
    double tyy;

    double t0xx;
    double t0xy;
    double t0yy;

    double tau;
    double tau_tau;
    double p_tr;
    double p_t;
    double S0;

    double dsjxx;
    double dsjxy;
    double dsjyy;

    double Ee_trxx;
    double Ee_trxy;
    double Ee_tryy;

    double Etauxx;
    double Etauxy;
    double Etauyy;

    double Np_trxx;
    double Np_trxy;
    double Np_tryy;

    double tau_star;

    double nup_tau;

    double dpxx;
    double dpxy;
    double dpyy;

    double const c = 1e-2;
    int i;

    for (i = 0; i < job->num_particles; i++) {
/*        fprintf(stderr, "%d: %f %f %f\n", i, job->particles[i].sxx,*/
/*            job->particles[i].sxy, job->particles[i].syy);*/
        if (job->particles[i].active == 0) {
            continue;
        }

        /* check if the density allows for supporting any stress */
        if ((job->particles[i].m / job->particles[i].v) < 1485.0f) {
            job->particles[i].sxx = 0;
            job->particles[i].sxy = 0;
            job->particles[i].syy = 0;
            continue;
        }

        /* calculate final total strain */
        Etauxx = job->particles[i].exx_t * job->dt + Exx;
        Etauxy = job->particles[i].exy_t * job->dt + Exy;
        Etauyy = job->particles[i].eyy_t * job->dt + Eyy;

        /* calculate trial elastic strain */
        Ee_trxx = Etauxx - Epxx;
        Ee_trxy = Etauxy - Epxy;
        Ee_tryy = Etauyy - Epyy;

        dsjxx = job->dt * (EMOD / (1 - NUMOD*NUMOD)) * ((job->particles[i].exx_t) + NUMOD * (job->particles[i].eyy_t));
        dsjxy = job->dt * (EMOD / (2 *(1 + NUMOD))) * (job->particles[i].exy_t);
        dsjyy = job->dt * (EMOD / (1 - NUMOD*NUMOD)) * ((job->particles[i].eyy_t) + NUMOD * (job->particles[i].exx_t));

        /* Jaumann spin terms */
        dsjxx -= 2 * job->dt * job->particles[i].wxy_t * job->particles[i].sxy;
        dsjxy += job->dt * job->particles[i].wxy_t * (job->particles[i].sxx - job->particles[i].syy);
        dsjyy += 2 * job->dt * job->particles[i].wxy_t * job->particles[i].sxy;

        /* calculate trial stress */
        txx = job->particles[i].sxx + dsjxx;
        txy = job->particles[i].sxy + dsjxy;
        tyy = job->particles[i].syy + dsjyy;

        /* pressure of trial  */
        p_tr = -0.5*(txx + tyy);

        if (p_tr < c) {
            job->particles[i].sxx = 0.5 * c / MU_S;
            job->particles[i].sxy = 0;
            job->particles[i].syy = 0.5 * c / MU_S;
            continue;
        }

        S0 = p_tr * MU_S;
        eta = B * GRAINS_D * sqrt(p_tr * GRAINS_RHO); 

        /* deviator of trial stress */
        t0xx = txx + p_tr;
        t0xy = txy;
        t0yy = tyy + p_tr;

        tau = sqrt(0.5f*(t0xx*t0xx + 2*t0xy*t0xy + t0yy*t0yy));

        if (tau > S0) {
            tau_tau = (tau + G * job->dt * S0 / eta) / (1.0 + G * job->dt / eta);
        } else {
            tau_tau = tau;
        }

        /* direction of plastic flow */
        Np_trxx = sqrt(0.5) * t0xx / tau;
        Np_trxy = sqrt(0.5) * t0xy / tau;
        Np_tryy = sqrt(0.5) * t0yy / tau;

        tau_star = tau_tau - S0;
        if (tau_star > 0) {
            nup_tau = tau_star / eta;
        } else {
            nup_tau = 0;
        }

        gammap += job->dt * nup_tau;

        /* plastic rate */
        dpxx = nup_tau * Np_trxx;
        dpxy = nup_tau * Np_trxy;
        dpyy = nup_tau * Np_tryy;

        Epxx += dpxx * job->dt;
        Epxy += dpxy * job->dt;
        Epyy += dpyy * job->dt;

/*        job->particles[i].sxx = txx - sqrt(2.0) * (tau - tau_tau) * Np_trxx;*/
/*        job->particles[i].sxy = txy - sqrt(2.0) * (tau - tau_tau) * Np_trxy;*/
/*        job->particles[i].syy = tyy - sqrt(2.0) * (tau - tau_tau) * Np_tryy;*/

        p_t = -0.5 * (job->particles[i].sxx + job->particles[i].syy);
        t0xx = job->particles[i].sxx + p_t;
        t0xy = job->particles[i].sxy;
        t0yy = job->particles[i].syy + p_t;
        tau = sqrt(0.5f*(t0xx*t0xx + 2*t0xy*t0xy + t0yy*t0yy));
        Np_trxx = sqrt(0.5) * t0xx / tau;
        Np_trxy = sqrt(0.5) * t0xy / tau;
        Np_tryy = sqrt(0.5) * t0yy / tau;
        
        eta = B * GRAINS_D * sqrt (GRAINS_RHO * p_t);
        if (tau > MU_S * p_t && p_t > c) {
            dpxx = sqrt(0.5) * ((tau - MU_S * p_t) / eta) * Np_trxx;
            dpxy = sqrt(0.5) * ((tau - MU_S * p_t) / eta) * Np_trxy;
            dpyy = sqrt(0.5) * ((tau - MU_S * p_t) / eta) * Np_tryy;
/*            fprintf(stderr, "dp != 0\n");*/
/*            dpxx = 0;*/
/*            dpxy = 0;*/
/*            dpyy = 0;*/
        } else {
            dpxx = 0;
            dpxy = 0;
            dpyy = 0;
        }

        if (job->t == 0 /* || job->step_number == 0 */) {
            dpxx = 0;
            dpxy = 0;
            dpyy = 0;
        }

        dsjxx = job->dt * (EMOD / (1 - NUMOD*NUMOD)) * ((job->particles[i].exx_t - dpxx) + NUMOD * (job->particles[i].eyy_t - dpyy));
        dsjxy = job->dt * (EMOD / (2 *(1 + NUMOD))) * (job->particles[i].exy_t - dpxy);
        dsjyy = job->dt * (EMOD / (1 - NUMOD*NUMOD)) * ((job->particles[i].eyy_t - dpyy) + NUMOD * (job->particles[i].exx_t - dpxx));
        dsjxx -= 2 * job->dt * job->particles[i].wxy_t * job->particles[i].sxy;
        dsjxy += job->dt * job->particles[i].wxy_t * (job->particles[i].sxx - job->particles[i].syy);
        dsjyy += 2 * job->dt * job->particles[i].wxy_t * job->particles[i].sxy;

        job->particles[i].sxx += dsjxx;
        job->particles[i].sxy += dsjxy;
        job->particles[i].syy += dsjyy;

        /* bulk viscosity */
    }

    return;
}
/*----------------------------------------------------------------------------*/

