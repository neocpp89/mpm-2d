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
#include "point.h"
#include "process.h"

#define G_MAG 1.0f
/*#define G_MAG 0.0f*/
#define RAMP_TIME 0.5f

/*----------------------------------------------------------------------------*/
void initial_loads(job_t *job)
{
    int i;

#ifdef __OEST__
    double y_max;

    for (i = 0; i < job->num_particles; i++) {
        job->particles[i].bx = 0;
        job->particles[i].by = -G_MAG;
    }

    y_max = 0;
    for (i = 0; i < job->num_particles; i++) {
        if (job->particles[i].y > y_max) {
            y_max = job->particles[i].y;
        }
    }

    for (i = 0; i < job->num_particles; i++) {
        job->particles[i].material = M_DRUCKER_PRAGER;
/*        printf("ymax = %lf\n", y_max);*/
        if (y_max - job->particles[i].y <= 1e-6) {
            job->particles[i].material = M_RIGID;
/*            job->particles[i].by = -10.0;*/
/*            printf("particle %d set as rigid. y = %lf\n",*/
/*                i, 1e-6);*/
        }
    }
#else
    for (i = 0; i < job->num_particles; i++) {
        job->particles[i].bx = 0;
        job->particles[i].by = 0;
        job->particles[i].material = M_DRUCKER_PRAGER;
    }
#endif

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void time_varying_loads(job_t *job)
{
    int i;
    double gravity;

    if ((RAMP_TIME - job->t) > 0) {
        gravity = -G_MAG * (job->t / RAMP_TIME);
    } else {
        gravity = -G_MAG;
    }
/*    printf("gravity set to %f.\n", gravity);*/
    for (i = 0; i < job->num_particles; i++) {
#ifdef __OEST__
        if (job->particles[i].material == M_RIGID) {
            job->particles[i].bx = 0;
            job->particles[i].by = gravity;
        } else {
            job->particles[i].bx = 0;
            job->particles[i].by = 0;
        }
#else
        job->particles[i].bx = 0;
        job->particles[i].by = gravity;
#endif
    }

    return;
}
/*----------------------------------------------------------------------------*/

