/**
    \file implicit.c
    \author Sachith Dunatunga
    \date 20.10.2014

    mpm_2d -- An implementation of the Material Point Method in 2D.
*/
#include <math.h>
#include "element.h"
#include "implicit.h"
#include "particle.h"
#include "process.h"
#include "material.h"

#define CHECK_ACTIVE(j,i) if (j->active[i] == 0) { continue; }

/*--build_elemental_stiffness-------------------------------------------------*/
void build_elemental_stiffness(job_t *job)
{
    /* clear all filled elements */
    for (size_t i = 0; i < job->num_elements; i++) {
        if (job->elements[i].filled == 0) {
            continue;
        }

        for (size_t j = 0; j < (NODAL_DOF * NODES_PER_ELEMENT); j++) {
            for (size_t k = 0; k < (NODAL_DOF * NODES_PER_ELEMENT); k++) {
                job->elements[i].kku_element[j][k] = 0;
            }
            job->elements[i].f_element[j] = 0;
        }
    }

    for (size_t i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);
        add_particle_stiffness(i, job);
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*---add_particle_stiffness---------------------------------------------------*/
void add_particle_stiffness(int idx, job_t *job)
{
    size_t i = idx;
    size_t p = job->in_element[i];

    /* --- F_int --- */
    job->elements[p].f_element[NODAL_DOF * 0 + 0] += ( job->particles[i].v*(job->b11[i]*job->particles[i].sxx + job->b21[i]*job->particles[i].sxy) );
    job->elements[p].f_element[NODAL_DOF * 0 + 1] += ( job->particles[i].v*(job->b11[i]*job->particles[i].sxy + job->b21[i]*job->particles[i].syy) );
    job->elements[p].f_element[NODAL_DOF * 1 + 0] += ( job->particles[i].v*(job->b12[i]*job->particles[i].sxx + job->b22[i]*job->particles[i].sxy) );
    job->elements[p].f_element[NODAL_DOF * 1 + 1] += ( job->particles[i].v*(job->b12[i]*job->particles[i].sxy + job->b22[i]*job->particles[i].syy) );
    job->elements[p].f_element[NODAL_DOF * 2 + 0] += ( job->particles[i].v*(job->b13[i]*job->particles[i].sxx + job->b23[i]*job->particles[i].sxy) );
    job->elements[p].f_element[NODAL_DOF * 2 + 1] += ( job->particles[i].v*(job->b13[i]*job->particles[i].sxy + job->b23[i]*job->particles[i].syy) );
    job->elements[p].f_element[NODAL_DOF * 3 + 0] += ( job->particles[i].v*(job->b14[i]*job->particles[i].sxx + job->b24[i]*job->particles[i].sxy) );
    job->elements[p].f_element[NODAL_DOF * 3 + 1] += ( job->particles[i].v*(job->b14[i]*job->particles[i].sxy + job->b24[i]*job->particles[i].syy) );

    /* --- K_mat --- Symmetric? True --- */
    job->elements[p].kku_element[NODAL_DOF * 0 + 0][NODAL_DOF * 0 + 0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*pow(job->b21[i], 2) - 1.0*pow(job->b11[i], 2) - 0.5*pow(job->b21[i], 2))/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 0 + 0][NODAL_DOF * 0 + 1] += ( 0.5*EMOD*job->b11[i]*job->b21[i]*job->particles[i].v/(-NUMOD + 1) );
    job->elements[p].kku_element[NODAL_DOF * 0 + 0][NODAL_DOF * 1 + 0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b21[i]*job->b22[i] - 1.0*job->b11[i]*job->b12[i] - 0.5*job->b21[i]*job->b22[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 0 + 0][NODAL_DOF * 1 + 1] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b11[i]*job->b22[i] + 0.5*NUMOD*job->b12[i]*job->b21[i] - 0.5*job->b12[i]*job->b21[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 0 + 0][NODAL_DOF * 2 + 0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b21[i]*job->b23[i] - 1.0*job->b11[i]*job->b13[i] - 0.5*job->b21[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 0 + 0][NODAL_DOF * 2 + 1] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b11[i]*job->b23[i] + 0.5*NUMOD*job->b13[i]*job->b21[i] - 0.5*job->b13[i]*job->b21[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 0 + 0][NODAL_DOF * 3 + 0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b21[i]*job->b24[i] - 1.0*job->b11[i]*job->b14[i] - 0.5*job->b21[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 0 + 0][NODAL_DOF * 3 + 1] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b11[i]*job->b24[i] + 0.5*NUMOD*job->b14[i]*job->b21[i] - 0.5*job->b14[i]*job->b21[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 0 + 1][NODAL_DOF * 0 + 0] += ( 0.5*EMOD*job->b11[i]*job->b21[i]*job->particles[i].v/(-NUMOD + 1) );
    job->elements[p].kku_element[NODAL_DOF * 0 + 1][NODAL_DOF * 0 + 1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*pow(job->b11[i], 2) - 0.5*pow(job->b11[i], 2) - 1.0*pow(job->b21[i], 2))/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 0 + 1][NODAL_DOF * 1 + 0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b22[i] - 1.0*NUMOD*job->b12[i]*job->b21[i] - 0.5*job->b11[i]*job->b22[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 0 + 1][NODAL_DOF * 1 + 1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b12[i] - 0.5*job->b11[i]*job->b12[i] - 1.0*job->b21[i]*job->b22[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 0 + 1][NODAL_DOF * 2 + 0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b23[i] - 1.0*NUMOD*job->b13[i]*job->b21[i] - 0.5*job->b11[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 0 + 1][NODAL_DOF * 2 + 1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b13[i] - 0.5*job->b11[i]*job->b13[i] - 1.0*job->b21[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 0 + 1][NODAL_DOF * 3 + 0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b24[i] - 1.0*NUMOD*job->b14[i]*job->b21[i] - 0.5*job->b11[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 0 + 1][NODAL_DOF * 3 + 1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b14[i] - 0.5*job->b11[i]*job->b14[i] - 1.0*job->b21[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 1 + 0][NODAL_DOF * 0 + 0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b21[i]*job->b22[i] - 1.0*job->b11[i]*job->b12[i] - 0.5*job->b21[i]*job->b22[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 1 + 0][NODAL_DOF * 0 + 1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b22[i] - 1.0*NUMOD*job->b12[i]*job->b21[i] - 0.5*job->b11[i]*job->b22[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 1 + 0][NODAL_DOF * 1 + 0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*pow(job->b22[i], 2) - 1.0*pow(job->b12[i], 2) - 0.5*pow(job->b22[i], 2))/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 1 + 0][NODAL_DOF * 1 + 1] += ( 0.5*EMOD*job->b12[i]*job->b22[i]*job->particles[i].v/(-NUMOD + 1) );
    job->elements[p].kku_element[NODAL_DOF * 1 + 0][NODAL_DOF * 2 + 0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b22[i]*job->b23[i] - 1.0*job->b12[i]*job->b13[i] - 0.5*job->b22[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 1 + 0][NODAL_DOF * 2 + 1] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b12[i]*job->b23[i] + 0.5*NUMOD*job->b13[i]*job->b22[i] - 0.5*job->b13[i]*job->b22[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 1 + 0][NODAL_DOF * 3 + 0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b22[i]*job->b24[i] - 1.0*job->b12[i]*job->b14[i] - 0.5*job->b22[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 1 + 0][NODAL_DOF * 3 + 1] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b12[i]*job->b24[i] + 0.5*NUMOD*job->b14[i]*job->b22[i] - 0.5*job->b14[i]*job->b22[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 1 + 1][NODAL_DOF * 0 + 0] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b11[i]*job->b22[i] + 0.5*NUMOD*job->b12[i]*job->b21[i] - 0.5*job->b12[i]*job->b21[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 1 + 1][NODAL_DOF * 0 + 1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b12[i] - 0.5*job->b11[i]*job->b12[i] - 1.0*job->b21[i]*job->b22[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 1 + 1][NODAL_DOF * 1 + 0] += ( 0.5*EMOD*job->b12[i]*job->b22[i]*job->particles[i].v/(-NUMOD + 1) );
    job->elements[p].kku_element[NODAL_DOF * 1 + 1][NODAL_DOF * 1 + 1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*pow(job->b12[i], 2) - 0.5*pow(job->b12[i], 2) - 1.0*pow(job->b22[i], 2))/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 1 + 1][NODAL_DOF * 2 + 0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b12[i]*job->b23[i] - 1.0*NUMOD*job->b13[i]*job->b22[i] - 0.5*job->b12[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 1 + 1][NODAL_DOF * 2 + 1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b12[i]*job->b13[i] - 0.5*job->b12[i]*job->b13[i] - 1.0*job->b22[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 1 + 1][NODAL_DOF * 3 + 0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b12[i]*job->b24[i] - 1.0*NUMOD*job->b14[i]*job->b22[i] - 0.5*job->b12[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 1 + 1][NODAL_DOF * 3 + 1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b12[i]*job->b14[i] - 0.5*job->b12[i]*job->b14[i] - 1.0*job->b22[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 2 + 0][NODAL_DOF * 0 + 0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b21[i]*job->b23[i] - 1.0*job->b11[i]*job->b13[i] - 0.5*job->b21[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 2 + 0][NODAL_DOF * 0 + 1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b23[i] - 1.0*NUMOD*job->b13[i]*job->b21[i] - 0.5*job->b11[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 2 + 0][NODAL_DOF * 1 + 0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b22[i]*job->b23[i] - 1.0*job->b12[i]*job->b13[i] - 0.5*job->b22[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 2 + 0][NODAL_DOF * 1 + 1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b12[i]*job->b23[i] - 1.0*NUMOD*job->b13[i]*job->b22[i] - 0.5*job->b12[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 2 + 0][NODAL_DOF * 2 + 0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*pow(job->b23[i], 2) - 1.0*pow(job->b13[i], 2) - 0.5*pow(job->b23[i], 2))/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 2 + 0][NODAL_DOF * 2 + 1] += ( 0.5*EMOD*job->b13[i]*job->b23[i]*job->particles[i].v/(-NUMOD + 1) );
    job->elements[p].kku_element[NODAL_DOF * 2 + 0][NODAL_DOF * 3 + 0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b23[i]*job->b24[i] - 1.0*job->b13[i]*job->b14[i] - 0.5*job->b23[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 2 + 0][NODAL_DOF * 3 + 1] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b13[i]*job->b24[i] + 0.5*NUMOD*job->b14[i]*job->b23[i] - 0.5*job->b14[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 2 + 1][NODAL_DOF * 0 + 0] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b11[i]*job->b23[i] + 0.5*NUMOD*job->b13[i]*job->b21[i] - 0.5*job->b13[i]*job->b21[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 2 + 1][NODAL_DOF * 0 + 1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b13[i] - 0.5*job->b11[i]*job->b13[i] - 1.0*job->b21[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 2 + 1][NODAL_DOF * 1 + 0] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b12[i]*job->b23[i] + 0.5*NUMOD*job->b13[i]*job->b22[i] - 0.5*job->b13[i]*job->b22[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 2 + 1][NODAL_DOF * 1 + 1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b12[i]*job->b13[i] - 0.5*job->b12[i]*job->b13[i] - 1.0*job->b22[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 2 + 1][NODAL_DOF * 2 + 0] += ( 0.5*EMOD*job->b13[i]*job->b23[i]*job->particles[i].v/(-NUMOD + 1) );
    job->elements[p].kku_element[NODAL_DOF * 2 + 1][NODAL_DOF * 2 + 1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*pow(job->b13[i], 2) - 0.5*pow(job->b13[i], 2) - 1.0*pow(job->b23[i], 2))/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 2 + 1][NODAL_DOF * 3 + 0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b13[i]*job->b24[i] - 1.0*NUMOD*job->b14[i]*job->b23[i] - 0.5*job->b13[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 2 + 1][NODAL_DOF * 3 + 1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b13[i]*job->b14[i] - 0.5*job->b13[i]*job->b14[i] - 1.0*job->b23[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 3 + 0][NODAL_DOF * 0 + 0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b21[i]*job->b24[i] - 1.0*job->b11[i]*job->b14[i] - 0.5*job->b21[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 3 + 0][NODAL_DOF * 0 + 1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b24[i] - 1.0*NUMOD*job->b14[i]*job->b21[i] - 0.5*job->b11[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 3 + 0][NODAL_DOF * 1 + 0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b22[i]*job->b24[i] - 1.0*job->b12[i]*job->b14[i] - 0.5*job->b22[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 3 + 0][NODAL_DOF * 1 + 1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b12[i]*job->b24[i] - 1.0*NUMOD*job->b14[i]*job->b22[i] - 0.5*job->b12[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 3 + 0][NODAL_DOF * 2 + 0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b23[i]*job->b24[i] - 1.0*job->b13[i]*job->b14[i] - 0.5*job->b23[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 3 + 0][NODAL_DOF * 2 + 1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b13[i]*job->b24[i] - 1.0*NUMOD*job->b14[i]*job->b23[i] - 0.5*job->b13[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 3 + 0][NODAL_DOF * 3 + 0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*pow(job->b24[i], 2) - 1.0*pow(job->b14[i], 2) - 0.5*pow(job->b24[i], 2))/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 3 + 0][NODAL_DOF * 3 + 1] += ( 0.5*EMOD*job->b14[i]*job->b24[i]*job->particles[i].v/(-NUMOD + 1) );
    job->elements[p].kku_element[NODAL_DOF * 3 + 1][NODAL_DOF * 0 + 0] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b11[i]*job->b24[i] + 0.5*NUMOD*job->b14[i]*job->b21[i] - 0.5*job->b14[i]*job->b21[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 3 + 1][NODAL_DOF * 0 + 1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b14[i] - 0.5*job->b11[i]*job->b14[i] - 1.0*job->b21[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 3 + 1][NODAL_DOF * 1 + 0] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b12[i]*job->b24[i] + 0.5*NUMOD*job->b14[i]*job->b22[i] - 0.5*job->b14[i]*job->b22[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 3 + 1][NODAL_DOF * 1 + 1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b12[i]*job->b14[i] - 0.5*job->b12[i]*job->b14[i] - 1.0*job->b22[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 3 + 1][NODAL_DOF * 2 + 0] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b13[i]*job->b24[i] + 0.5*NUMOD*job->b14[i]*job->b23[i] - 0.5*job->b14[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 3 + 1][NODAL_DOF * 2 + 1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b13[i]*job->b14[i] - 0.5*job->b13[i]*job->b14[i] - 1.0*job->b23[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[NODAL_DOF * 3 + 1][NODAL_DOF * 3 + 0] += ( 0.5*EMOD*job->b14[i]*job->b24[i]*job->particles[i].v/(-NUMOD + 1) );
    job->elements[p].kku_element[NODAL_DOF * 3 + 1][NODAL_DOF * 3 + 1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*pow(job->b14[i], 2) - 0.5*pow(job->b14[i], 2) - 1.0*pow(job->b24[i], 2))/(pow(NUMOD, 2) - 1) );

    /* --- K_geo --- Symmetric? True --- */
    job->elements[p].kku_element[NODAL_DOF * 0 + 0][NODAL_DOF * 0 + 0] += ( job->particles[i].v*(pow(job->b11[i], 2)*job->particles[i].sxx + 2*job->b11[i]*job->b21[i]*job->particles[i].sxy + pow(job->b21[i], 2)*job->particles[i].syy) );
    job->elements[p].kku_element[NODAL_DOF * 0 + 0][NODAL_DOF * 1 + 0] += ( job->particles[i].v*(job->b11[i]*job->b12[i]*job->particles[i].sxx + job->b11[i]*job->b22[i]*job->particles[i].sxy + job->b12[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b22[i]*job->particles[i].syy) );
    job->elements[p].kku_element[NODAL_DOF * 0 + 0][NODAL_DOF * 2 + 0] += ( job->particles[i].v*(job->b11[i]*job->b13[i]*job->particles[i].sxx + job->b11[i]*job->b23[i]*job->particles[i].sxy + job->b13[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b23[i]*job->particles[i].syy) );
    job->elements[p].kku_element[NODAL_DOF * 0 + 0][NODAL_DOF * 3 + 0] += ( job->particles[i].v*(job->b11[i]*job->b14[i]*job->particles[i].sxx + job->b11[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b24[i]*job->particles[i].syy) );
    job->elements[p].kku_element[NODAL_DOF * 0 + 1][NODAL_DOF * 0 + 1] += ( job->particles[i].v*(pow(job->b11[i], 2)*job->particles[i].sxx + 2*job->b11[i]*job->b21[i]*job->particles[i].sxy + pow(job->b21[i], 2)*job->particles[i].syy) );
    job->elements[p].kku_element[NODAL_DOF * 0 + 1][NODAL_DOF * 1 + 1] += ( job->particles[i].v*(job->b11[i]*job->b12[i]*job->particles[i].sxx + job->b11[i]*job->b22[i]*job->particles[i].sxy + job->b12[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b22[i]*job->particles[i].syy) );
    job->elements[p].kku_element[NODAL_DOF * 0 + 1][NODAL_DOF * 2 + 1] += ( job->particles[i].v*(job->b11[i]*job->b13[i]*job->particles[i].sxx + job->b11[i]*job->b23[i]*job->particles[i].sxy + job->b13[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b23[i]*job->particles[i].syy) );
    job->elements[p].kku_element[NODAL_DOF * 0 + 1][NODAL_DOF * 3 + 1] += ( job->particles[i].v*(job->b11[i]*job->b14[i]*job->particles[i].sxx + job->b11[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b24[i]*job->particles[i].syy) );
    job->elements[p].kku_element[NODAL_DOF * 1 + 0][NODAL_DOF * 0 + 0] += ( job->particles[i].v*(job->b11[i]*job->b12[i]*job->particles[i].sxx + job->b11[i]*job->b22[i]*job->particles[i].sxy + job->b12[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b22[i]*job->particles[i].syy) );
    job->elements[p].kku_element[NODAL_DOF * 1 + 0][NODAL_DOF * 1 + 0] += ( job->particles[i].v*(pow(job->b12[i], 2)*job->particles[i].sxx + 2*job->b12[i]*job->b22[i]*job->particles[i].sxy + pow(job->b22[i], 2)*job->particles[i].syy) );
    job->elements[p].kku_element[NODAL_DOF * 1 + 0][NODAL_DOF * 2 + 0] += ( job->particles[i].v*(job->b12[i]*job->b13[i]*job->particles[i].sxx + job->b12[i]*job->b23[i]*job->particles[i].sxy + job->b13[i]*job->b22[i]*job->particles[i].sxy + job->b22[i]*job->b23[i]*job->particles[i].syy) );
    job->elements[p].kku_element[NODAL_DOF * 1 + 0][NODAL_DOF * 3 + 0] += ( job->particles[i].v*(job->b12[i]*job->b14[i]*job->particles[i].sxx + job->b12[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b22[i]*job->particles[i].sxy + job->b22[i]*job->b24[i]*job->particles[i].syy) );
    job->elements[p].kku_element[NODAL_DOF * 1 + 1][NODAL_DOF * 0 + 1] += ( job->particles[i].v*(job->b11[i]*job->b12[i]*job->particles[i].sxx + job->b11[i]*job->b22[i]*job->particles[i].sxy + job->b12[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b22[i]*job->particles[i].syy) );
    job->elements[p].kku_element[NODAL_DOF * 1 + 1][NODAL_DOF * 1 + 1] += ( job->particles[i].v*(pow(job->b12[i], 2)*job->particles[i].sxx + 2*job->b12[i]*job->b22[i]*job->particles[i].sxy + pow(job->b22[i], 2)*job->particles[i].syy) );
    job->elements[p].kku_element[NODAL_DOF * 1 + 1][NODAL_DOF * 2 + 1] += ( job->particles[i].v*(job->b12[i]*job->b13[i]*job->particles[i].sxx + job->b12[i]*job->b23[i]*job->particles[i].sxy + job->b13[i]*job->b22[i]*job->particles[i].sxy + job->b22[i]*job->b23[i]*job->particles[i].syy) );
    job->elements[p].kku_element[NODAL_DOF * 1 + 1][NODAL_DOF * 3 + 1] += ( job->particles[i].v*(job->b12[i]*job->b14[i]*job->particles[i].sxx + job->b12[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b22[i]*job->particles[i].sxy + job->b22[i]*job->b24[i]*job->particles[i].syy) );
    job->elements[p].kku_element[NODAL_DOF * 2 + 0][NODAL_DOF * 0 + 0] += ( job->particles[i].v*(job->b11[i]*job->b13[i]*job->particles[i].sxx + job->b11[i]*job->b23[i]*job->particles[i].sxy + job->b13[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b23[i]*job->particles[i].syy) );
    job->elements[p].kku_element[NODAL_DOF * 2 + 0][NODAL_DOF * 1 + 0] += ( job->particles[i].v*(job->b12[i]*job->b13[i]*job->particles[i].sxx + job->b12[i]*job->b23[i]*job->particles[i].sxy + job->b13[i]*job->b22[i]*job->particles[i].sxy + job->b22[i]*job->b23[i]*job->particles[i].syy) );
    job->elements[p].kku_element[NODAL_DOF * 2 + 0][NODAL_DOF * 2 + 0] += ( job->particles[i].v*(pow(job->b13[i], 2)*job->particles[i].sxx + 2*job->b13[i]*job->b23[i]*job->particles[i].sxy + pow(job->b23[i], 2)*job->particles[i].syy) );
    job->elements[p].kku_element[NODAL_DOF * 2 + 0][NODAL_DOF * 3 + 0] += ( job->particles[i].v*(job->b13[i]*job->b14[i]*job->particles[i].sxx + job->b13[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b23[i]*job->particles[i].sxy + job->b23[i]*job->b24[i]*job->particles[i].syy) );
    job->elements[p].kku_element[NODAL_DOF * 2 + 1][NODAL_DOF * 0 + 1] += ( job->particles[i].v*(job->b11[i]*job->b13[i]*job->particles[i].sxx + job->b11[i]*job->b23[i]*job->particles[i].sxy + job->b13[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b23[i]*job->particles[i].syy) );
    job->elements[p].kku_element[NODAL_DOF * 2 + 1][NODAL_DOF * 1 + 1] += ( job->particles[i].v*(job->b12[i]*job->b13[i]*job->particles[i].sxx + job->b12[i]*job->b23[i]*job->particles[i].sxy + job->b13[i]*job->b22[i]*job->particles[i].sxy + job->b22[i]*job->b23[i]*job->particles[i].syy) );
    job->elements[p].kku_element[NODAL_DOF * 2 + 1][NODAL_DOF * 2 + 1] += ( job->particles[i].v*(pow(job->b13[i], 2)*job->particles[i].sxx + 2*job->b13[i]*job->b23[i]*job->particles[i].sxy + pow(job->b23[i], 2)*job->particles[i].syy) );
    job->elements[p].kku_element[NODAL_DOF * 2 + 1][NODAL_DOF * 3 + 1] += ( job->particles[i].v*(job->b13[i]*job->b14[i]*job->particles[i].sxx + job->b13[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b23[i]*job->particles[i].sxy + job->b23[i]*job->b24[i]*job->particles[i].syy) );
    job->elements[p].kku_element[NODAL_DOF * 3 + 0][NODAL_DOF * 0 + 0] += ( job->particles[i].v*(job->b11[i]*job->b14[i]*job->particles[i].sxx + job->b11[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b24[i]*job->particles[i].syy) );
    job->elements[p].kku_element[NODAL_DOF * 3 + 0][NODAL_DOF * 1 + 0] += ( job->particles[i].v*(job->b12[i]*job->b14[i]*job->particles[i].sxx + job->b12[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b22[i]*job->particles[i].sxy + job->b22[i]*job->b24[i]*job->particles[i].syy) );
    job->elements[p].kku_element[NODAL_DOF * 3 + 0][NODAL_DOF * 2 + 0] += ( job->particles[i].v*(job->b13[i]*job->b14[i]*job->particles[i].sxx + job->b13[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b23[i]*job->particles[i].sxy + job->b23[i]*job->b24[i]*job->particles[i].syy) );
    job->elements[p].kku_element[NODAL_DOF * 3 + 0][NODAL_DOF * 3 + 0] += ( job->particles[i].v*(pow(job->b14[i], 2)*job->particles[i].sxx + 2*job->b14[i]*job->b24[i]*job->particles[i].sxy + pow(job->b24[i], 2)*job->particles[i].syy) );
    job->elements[p].kku_element[NODAL_DOF * 3 + 1][NODAL_DOF * 0 + 1] += ( job->particles[i].v*(job->b11[i]*job->b14[i]*job->particles[i].sxx + job->b11[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b24[i]*job->particles[i].syy) );
    job->elements[p].kku_element[NODAL_DOF * 3 + 1][NODAL_DOF * 1 + 1] += ( job->particles[i].v*(job->b12[i]*job->b14[i]*job->particles[i].sxx + job->b12[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b22[i]*job->particles[i].sxy + job->b22[i]*job->b24[i]*job->particles[i].syy) );
    job->elements[p].kku_element[NODAL_DOF * 3 + 1][NODAL_DOF * 2 + 1] += ( job->particles[i].v*(job->b13[i]*job->b14[i]*job->particles[i].sxx + job->b13[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b23[i]*job->particles[i].sxy + job->b23[i]*job->b24[i]*job->particles[i].syy) );
    job->elements[p].kku_element[NODAL_DOF * 3 + 1][NODAL_DOF * 3 + 1] += ( job->particles[i].v*(pow(job->b14[i], 2)*job->particles[i].sxx + 2*job->b14[i]*job->b24[i]*job->particles[i].sxy + pow(job->b24[i], 2)*job->particles[i].syy) );

    return;
}
/*----------------------------------------------------------------------------*/

