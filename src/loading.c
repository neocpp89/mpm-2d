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

#define G_MAG 9.80

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
    size_t i, j; 
    double gravity;
    double compress;
    double p;

    int sqrtnp;
    int sqrtnel;

#if 0
    p = 0;
    for (i = 0; i < job->num_particles; i++) {
        p += -0.5*(job->particles[i].sxx + job->particles[i].syy) / job->num_particles;
    }
/*    printf("[0].pressure = %f\n", p);*/

    /* we compressed it to about 400 Pa */
    if (job->step_number == 0 && p >= 100) {
        job->step_number++;
        job->step_start_time = job->t;
        fprintf(job->output.log_fd, "Done compressing laterally. Moving to next step.\n");
        fflush(job->output.log_fd);
    }

    if (job->step_number == 1 && (job->t >= job->step_start_time + 0.0)) {
        job->step_number++;
        job->step_start_time = job->t;
        fprintf(job->output.log_fd, "Done waiting. Moving to next step.\n");
        fflush(job->output.log_fd);
        /* deactivate particles in 'frozen' elements (left vertical column). */
        sqrtnel = (int)sqrt(job->num_elements);
        for (i = 0; i < job->num_particles; i++) {
/*            for (j = 0; j < sqrtnel; j++) {*/
                if (job->in_element[i] % sqrtnel == 0) {
                    job->active[i] = 0;
                    fprintf(job->output.log_fd,
                        "Particle %d @(%g, %g) in frozen element, deactivating.\n",
                        i, job->particles[i].x, job->particles[i].y);
                }
/*            }*/
        }
    }
#endif

    switch (job->step_number) {
#if 0
        case 0:
            sqrtnp = (int)sqrt(job->num_particles);
            compress = 100 * (job->t - job->step_start_time);
            if (compress > 100) { compress = 100; }
            
            for (i = 0; i < sqrtnp; i++) {
                job->particles[i].bx = compress;
                job->particles[i].by = 0;
            }
        break;
        case 1:
            /* just wait a bit for waves to die out. */
            for (i = 0; i < job->num_particles; i++) {
                job->particles[i].bx = 0;
                job->particles[i].by = 0;
            }
        break;
        case 2:
#endif
        default:
            gravity = -G_MAG;

            for (i = 0; i < job->num_particles; i++) {
                job->particles[i].bx = 0;
                job->particles[i].by = gravity;
            }
        break;
    }

    return;
}
/*----------------------------------------------------------------------------*/

