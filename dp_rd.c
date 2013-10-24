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

#define EMOD 1e6
#define NUMOD 0.3

#define G (EMOD / (2.0f * (1.0f + NUMOD)))
#define K (EMOD / (3.0f * (1.0f - 2*NUMOD)))

#define PI 3.1415926535897932384626433
#define PHI (30.0*PI/180.0)

#define MU_S 0.577350269
#define GRAINS_RHO 3000
#define GRAINS_D 0.01

#define Epxx jp(state[0])
#define Epxy jp(state[1])
#define Epyy jp(state[2])
#define gf jp(state[3])
#define eta jp(state[4])
#define beta jp(state[5])
#define gammap jp(state[9])
#define Ef_mag jp(state[10])
#define Etxx jp(state[6])
#define Etxy jp(state[7])
#define Etyy jp(state[8])

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
        Etxx = 0;
        Etxy = 0;
        Etyy = 0;
        eta = 0;
        gf = 0;
    }

    fprintf(stderr, "%s:%s: done initializing material.\n", __FILE__,  __func__);
    return;
}
/*----------------------------------------------------------------------------*/

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

    double dp_xx;
    double dp_xy;
    double dp_yy;

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

/*        dsjxx -= 2 * job->dt * job->particles[i].wxy_t * job->particles[i].real_sxy;*/
/*        dsjxy += job->dt * job->particles[i].wxy_t * (job->particles[i].real_sxx - job->particles[i].real_syy);*/
/*        dsjyy += 2 * job->dt * job->particles[i].wxy_t * job->particles[i].real_sxy;*/

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
            continue;
        }

        if ((job->particles[i].m / job->particles[i].v) < 1485.0f) {
            density_flag = 1;
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
        } else {
            dp_xx = 1e-1 * f * t0xx / qm;
            dp_xy = 1e-1 * f * t0xy / qm;
            dp_yy = 1e-1 * f * t0yy / qm;
            job->particles[i].sxx = txx - job->dt * (EMOD / (1 - NUMOD*NUMOD)) * ((dp_xx) + NUMOD * (dp_yy));
            job->particles[i].sxy = txy - job->dt * (EMOD / (2 *(1 + NUMOD))) * (dp_xy);
            job->particles[i].syy = tyy - job->dt * (EMOD / (1 - NUMOD*NUMOD)) * ((dp_yy) + NUMOD * (dp_xx));
            Epxx = Epxx + dp_xx * job->dt;
            Epxy = Epxy + dp_xy * job->dt;
            Epyy = Epyy + dp_yy * job->dt;
            gammap = sqrt(Epxx*Epxx + 2*Epxy*Epxy + Epyy*Epyy);
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/
