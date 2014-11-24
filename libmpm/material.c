/**
    \file material.c
    \author Sachith Dunatunga
    \date 04.06.12

    mpm_2d -- Linear elastic stub for a material model in MPM.
*/
#include "particle.h"
#include "process.h"

// contains EMOD, NUMOD
#include "material.h"

// special state variable for equivalent plastic shear strain rate.
#define gammap state[9]

/*----------------------------------------------------------------------------*/
void material_init_linear_elastic(job_t *job)
{
    // clear state variables
    for (size_t i = 0; i < job->num_particles; i++) {
        for (size_t j = 0; j < DEPVAR; j++) {
            job->particles[i].state[j] = 0;
        }
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

    for (size_t i = 0; i < job->num_particles; i++) {
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

