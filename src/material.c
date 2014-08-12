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
void material_init_linear_elastic(job_t *job)
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
/*        Exx = 0;*/
/*        Exy = 0;*/
/*        Eyy = 0;*/
        gf = 0;
        gammap = 0;
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void calculate_stress_linear_elastic(job_t *job)
{
    double dsjxx;
    double dsjxy;
    double dsjyy;

    int i;

    for (i = 0; i < job->num_particles; i++) {
        if (job->active[i] == 0) {
            continue;
        }

        dsjxx = job->dt * (EMOD / (1 - NUMOD*NUMOD)) * ((job->particles[i].exx_t) + NUMOD * (job->particles[i].eyy_t));
        dsjxy = job->dt * (EMOD / (2 *(1 + NUMOD))) * (job->particles[i].exy_t);
        dsjyy = job->dt * (EMOD / (1 - NUMOD*NUMOD)) * ((job->particles[i].eyy_t) + NUMOD * (job->particles[i].exx_t));

        dsjxx += 2 * job->dt * job->particles[i].wxy_t * job->particles[i].sxy;
        dsjxy -= job->dt * job->particles[i].wxy_t * (job->particles[i].sxx - job->particles[i].syy);
        dsjyy -= 2 * job->dt * job->particles[i].wxy_t * job->particles[i].sxy;

        job->particles[i].sxx = job->particles[i].sxx + dsjxx;
        job->particles[i].sxy = job->particles[i].sxy + dsjxy;
        job->particles[i].syy = job->particles[i].syy + dsjyy;
    }

    return;
}
/*----------------------------------------------------------------------------*/


