/**
    \file loading.c
    \author Sachith Dunatunga
    \date 21.06.13

    mpm_2d -- An implementation of the Material Point Method in 2D.
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "material.h"
#include "particle.h"
#include "process.h"

#define G_MAG 9.80
#define RAMP_TIME 0.0

/*----------------------------------------------------------------------------*/
void initial_loads(job_t *job)
{
    for (size_t i = 0; i < job->num_particles; i++) {
        job->particles[i].bx = 0;
        job->particles[i].by = 0;
        job->particles[i].material = M_DRUCKER_PRAGER;
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void time_varying_loads(job_t *job)
{
    double gravity;

    switch (job->step_number) {
        default:
            if ((RAMP_TIME + job->step_start_time) > job->t) {
                gravity = -G_MAG * ((job->t - job->step_start_time) / RAMP_TIME);
            } else {
                gravity = -G_MAG;
            }

            for (size_t i = 0; i < job->num_particles; i++) {
                job->particles[i].bx = 0;
                job->particles[i].by = gravity;
            }
        break;
    }

    return;
}
/*----------------------------------------------------------------------------*/

