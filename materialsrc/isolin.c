/**
    \file isolin.c
    \author Sachith Dunatunga
    \date 23.10.13

    The isotropic linear elastic material model.
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

#define MAT_VERSION_STRING "1.0" __DATE__ " " __TIME__

/* Material Constants (set by configuration file). */
double E, nu, G, K;
double lambda;

void material_init(job_t *job);
void calculate_stress_threaded(threadtask_t *task);
void calculate_stress(job_t *job);

/*----------------------------------------------------------------------------*/
void material_init(job_t *job)
{
    size_t i, j;

    for (i = 0; i < job->num_particles; i++) {
        for (j = 0; j < DEPVAR; j++) {
            job->particles[i].state[j] = 0;
        }
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
        lambda = K - 2.0 * G / 3.0;
        printf("%s:%s: properties (E = %g, nu = %g, G = %g, K = %g).\n",
            __FILE__, __func__, E, nu, G , K);
    }

    printf("%s:%s: (material version %s) done initializing material.\n",
        __FILE__,  __func__, MAT_VERSION_STRING);
    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
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

/*----------------------------------------------------------------------------*/
void calculate_stress_threaded(threadtask_t *task)
{
    job_t *job = task->job;

    /* Since this is local, we can split the particles among the threads. */
    size_t p_start = task->offset;
    size_t p_stop = task->offset + task->blocksize;

    double dsjxx;
    double dsjxy;
    double dsjyy;
    double evoldot;
    double trD;
    size_t i;

    for (i = p_start; i < p_stop; i++) {
        if (job->active[i] == 0) {
            continue;
        }

        trD = job->particles[i].exx_t + job->particles[i].eyy_t;
        dsjxx = lambda * trD + 2.0 * G * job->particles[i].exx_t;
        dsjxy = 2.0 * G * job->particles[i].exy_t;
        dsjyy = lambda * trD + 2.0 * G * job->particles[i].eyy_t;

        dsjxx -= 2 * job->dt * job->particles[i].wxy_t * job->particles[i].sxy;
        dsjxy += job->dt * job->particles[i].wxy_t * (job->particles[i].sxx - job->particles[i].syy);
        dsjyy += 2 * job->dt * job->particles[i].wxy_t * job->particles[i].sxy;

/*
        evoldot = (job->particles[i].exx_t + job->particles[i].eyy_t);
        dsjxx += 1e-4 * K * evoldot;
        dsjyy += 1e-4 * K * evoldot;
*/

        job->particles[i].sxx += dsjxx;
        job->particles[i].sxy += dsjxy;
        job->particles[i].syy += dsjyy;
    }

    return;
}
/*----------------------------------------------------------------------------*/
