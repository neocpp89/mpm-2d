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
#include "interpolate.h"
#include "particle.h"
#include "node.h"
#include "process.h"
#include "material.h"

#define signum(x) ((int)((0 < x) - (x < 0)))

#define jp(x) job->particles[i].x
#define SVISC 1e3
#define SYIELD 1e2

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

#define EMOD 1e6
#define NUMOD 0.3

#define G (EMOD / (2.0f * (1.0f + NUMOD)))
#define K (EMOD / (3.0f * (1.0f - 2*NUMOD)))

#define PI 3.1415926535897932384626433
#define PHI (30.0*PI/180.0)

#define Epxx jp(state[0])
#define Epxy jp(state[1])
#define Epyy jp(state[2])
#define mu jp(state[3])
#define eta jp(state[4])
#define beta jp(state[5])
/*#define Exx jp(state[6])*/
/*#define Exy jp(state[7])*/
/*#define Eyy jp(state[8])*/
#define gammap jp(state[9])
#define Ef_mag jp(state[10])
#define gf jp(state[6])

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
        mu = mu0;
        eta = eta0;
        beta = hbeta*(mu0 - mucv);
        gf = 0;
        gammap = 0;
    }

    return;
}
/*----------------------------------------------------------------------------*/

/* Rate Independent Drucker-Prager (return along shear direction only). */
/*----------------------------------------------------------------------------*/
void calculate_stress(job_t *job)
{
    double f;
    double txx;
    double txy;
    double tyy;
    double t0xx;
    double t0xy;
    double t0yy;
    double qm;
    double pm;
    double m;
    double q_adj;

    double dsjxx;
    double dsjxy;
    double dsjyy;

    double density_flag;

    double const c = 1e-2;
    int i;

    for (i = 0; i < job->num_particles; i++) {
        if (job->particles[i].active == 0) {
            continue;
        }

        dsjxx = job->dt * (EMOD / (1 - NUMOD*NUMOD)) * ((job->particles[i].exx_t) + NUMOD * (job->particles[i].eyy_t));
        dsjxy = job->dt * (EMOD / (2 *(1 + NUMOD))) * (job->particles[i].exy_t);
        dsjyy = job->dt * (EMOD / (1 - NUMOD*NUMOD)) * ((job->particles[i].eyy_t) + NUMOD * (job->particles[i].exx_t));

        dsjxx -= 2 * job->dt * job->particles[i].wxy_t * job->particles[i].sxy;
        dsjxy += job->dt * job->particles[i].wxy_t * (job->particles[i].sxx - job->particles[i].syy);
        dsjyy += 2 * job->dt * job->particles[i].wxy_t * job->particles[i].sxy;

/*        dsjxx += 2 * job->dt * job->particles[i].wxy_t * job->particles[i].sxy;*/
/*        dsjxy -= job->dt * job->particles[i].wxy_t * (job->particles[i].sxx - job->particles[i].syy);*/
/*        dsjyy -= 2 * job->dt * job->particles[i].wxy_t * job->particles[i].sxy;*/

/*        dsjxx += 1e2 * (jp(exx_t) + jp(eyy_t));*/
/*        dsjyy += 1e2 * (jp(eyy_t) + jp(exx_t));*/

        txx = job->particles[i].sxx + dsjxx;
        txy = job->particles[i].sxy + dsjxy;
        tyy = job->particles[i].syy + dsjyy;
        pm = -0.5f*(txx + tyy);
        t0xx = txx + pm;
        t0xy = txy;
        t0yy = tyy + pm;
        qm = sqrt(0.5f*(t0xx*t0xx + 2*t0xy*t0xy + t0yy*t0yy));
        m = tan(PHI);

        f = qm - m*pm - c;

        if(job->particles[i].material == M_RIGID) {
            job->particles[i].sxx = txx;
            job->particles[i].sxy = txy;
            job->particles[i].syy = tyy;
/*            printf("particle %d is rigid.\n", i);*/
            continue;
        }

        if ((job->particles[i].m / job->particles[i].v) < 1485.0f) {
            density_flag = 1;
/*            printf("%4d: density %lf\n", i, (job->particles[i].m / job->particles[i].v));*/
        } else {
            density_flag = 0;
        }

        if (density_flag) {
            job->particles[i].sxx = 0;
            job->particles[i].sxy = 0;
            job->particles[i].syy = 0;
        } else if (f < 0) {
            job->particles[i].sxx = txx;
            job->particles[i].sxy = txy;
            job->particles[i].syy = tyy;
            job->particles[i].color = 1;
        } else if (f >= 0 && pm > -c/m) {
            Epxx += (f / qm) * t0xx / (2 * G);
            Epxy += (f / qm) * t0xy / (2 * G);
            Epyy += (f / qm) * t0yy / (2 * G);
            q_adj = m*pm + c;
            job->particles[i].sxx = (q_adj / qm) * t0xx - pm;
            job->particles[i].sxy = (q_adj / qm) * t0xy;
            job->particles[i].syy = (q_adj / qm) * t0yy - pm;
            job->particles[i].color = 2;
            gammap = sqrt(Epxx*Epxx + 2*Epxy*Epxy + Epyy*Epyy);
        } else if (pm <= -c/m) {
            job->particles[i].sxx = -0.5 * c / m;
            job->particles[i].sxy = 0;
            job->particles[i].syy = -0.5 * c / m;
            job->particles[i].color = 3;
        } else {
            fprintf(stderr, "u");
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/

