/**
    \file element.c
    \author Sachith Dunatunga
    \date 04.06.12

    mpm_2d -- An implementation of the Material Point Method in 2D.
*/
#include <math.h>
#include "element.h"
#include "particle.h"
#include "process.h"
#include "material.h"

#define CHECK_ACTIVE(j,i) if (j->particles[i].active == 0) { continue; }

void build_elemental_stiffness(job_t *job);
void add_particle_stiffness(int idx, job_t *job);


/*---node_number_to_coords----------------------------------------------------*/
void inline node_number_to_coords(double *x, double *y, int num, int N, double h)
{
    int i;
    int j;

    i = num % N;
    j = num / N;

    *x = i*h;
    *y = j*h;

    return;
}
/*----------------------------------------------------------------------------*/

/*--build_elemental_stiffness-------------------------------------------------*/
void build_elemental_stiffness(job_t *job)
{
    int i, j, k;

    /* clear all filled elements */
    for (i = 0; i < job->num_elements; i++) {
        if (job->elements[i].filled == 0) {
            continue;
        }

        for (j = 0; j < (NODAL_DOF * NODES_PER_ELEMENT); j++) {
            for (k = 0; k < (NODAL_DOF * NODES_PER_ELEMENT); k++) {
                job->elements[i].kku_element[j][k] = 0;
            }
            job->elements[i].f_element[j] = 0;
        }
    }

    for (i = 0; i < job->num_particles; i++) {
        CHECK_ACTIVE(job, i);
        add_particle_stiffness(i, job);
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*---add_particle_stiffness---------------------------------------------------*/
void add_particle_stiffness(int idx, job_t *job)
{
    int p;
    int i, j, k;
    double b[4][2];
    double h[4];
    double xi = 0.5;
    double g_loc;
    double gammadot;
    double tau;
    double pr;

    i = idx;

    b[0][0] = job->b11[i];
    b[1][0] = job->b12[i];
    b[2][0] = job->b13[i];
    b[3][0] = job->b14[i];
    b[0][1] = job->b21[i];
    b[1][1] = job->b22[i];
    b[2][1] = job->b23[i];
    b[3][1] = job->b24[i];
    h[0] = job->h1[i];
    h[1] = job->h2[i];
    h[2] = job->h3[i];
    h[3] = job->h4[i];

    p = job->in_element[i];

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


    /* nonlocal variable g (stored in particle.state[6]) */
    pr = -0.5 * (job->particles[i].sxx + job->particles[i].syy);
    tau = sqrt(0.5f*((job->particles[i].sxx + pr)*(job->particles[i].sxx + pr)
        + 2*job->particles[i].sxy*job->particles[i].sxy
        + (job->particles[i].syy + pr)*(job->particles[i].syy + pr)));
    gammadot = sqrt(2.0*(0.25*(job->particles[i].exx_t - job->particles[i].eyy_t)*(job->particles[i].exx_t - job->particles[i].eyy_t) +
        2*job->particles[i].exy_t*job->particles[i].exy_t +
        0.25*(job->particles[i].eyy_t - job->particles[i].exx_t)*(job->particles[i].eyy_t - job->particles[i].exx_t)
    ));
    g_loc = pr * gammadot / tau;

    for (j = 0; j < NODES_PER_ELEMENT; j++) {
        for (k = 0; k < NODES_PER_ELEMENT; k++) {
            job->elements[p].kku_element[NODAL_DOF * j + 2][NODAL_DOF * k + 2] +=
                job->particles[i].v * ((b[j][0]*b[k][0] + b[j][1]*b[k][1]) + h[j]*h[k]/(xi * xi));
        }
    }

    for (j = 0; j < NODES_PER_ELEMENT; j++) {
        job->elements[p].f_element[NODAL_DOF * j + 2] += 
            job->particles[i].v * (job->particles[i].state[6]*(b[j][0]*b[j][0] + b[j][1]*b[j][1]) + h[j]*((job->particles[i].state[6] - g_loc)/(xi * xi)));
    }

    return;
}
/*----------------------------------------------------------------------------*/

