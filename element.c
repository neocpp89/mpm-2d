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
    int i;

    i = idx;
    p = job->in_element[i];

    /* --- F_int --- */
    job->elements[p].f_element[0] += ( job->particles[i].v*(job->b11[i]*job->particles[i].sxx + job->b21[i]*job->particles[i].sxy) );
    job->elements[p].f_element[1] += ( job->particles[i].v*(job->b11[i]*job->particles[i].sxy + job->b21[i]*job->particles[i].syy) );
    job->elements[p].f_element[2] += ( job->particles[i].v*(job->b12[i]*job->particles[i].sxx + job->b22[i]*job->particles[i].sxy) );
    job->elements[p].f_element[3] += ( job->particles[i].v*(job->b12[i]*job->particles[i].sxy + job->b22[i]*job->particles[i].syy) );
    job->elements[p].f_element[4] += ( job->particles[i].v*(job->b13[i]*job->particles[i].sxx + job->b23[i]*job->particles[i].sxy) );
    job->elements[p].f_element[5] += ( job->particles[i].v*(job->b13[i]*job->particles[i].sxy + job->b23[i]*job->particles[i].syy) );
    job->elements[p].f_element[6] += ( job->particles[i].v*(job->b14[i]*job->particles[i].sxx + job->b24[i]*job->particles[i].sxy) );
    job->elements[p].f_element[7] += ( job->particles[i].v*(job->b14[i]*job->particles[i].sxy + job->b24[i]*job->particles[i].syy) );

    /* --- K_mat --- Symmetric? True --- */
    job->elements[p].kku_element[0][0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*pow(job->b21[i], 2) - 1.0*pow(job->b11[i], 2) - 0.5*pow(job->b21[i], 2))/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[0][1] += ( 0.5*EMOD*job->b11[i]*job->b21[i]*job->particles[i].v/(-NUMOD + 1) );
    job->elements[p].kku_element[0][2] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b21[i]*job->b22[i] - 1.0*job->b11[i]*job->b12[i] - 0.5*job->b21[i]*job->b22[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[0][3] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b11[i]*job->b22[i] + 0.5*NUMOD*job->b12[i]*job->b21[i] - 0.5*job->b12[i]*job->b21[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[0][4] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b21[i]*job->b23[i] - 1.0*job->b11[i]*job->b13[i] - 0.5*job->b21[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[0][5] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b11[i]*job->b23[i] + 0.5*NUMOD*job->b13[i]*job->b21[i] - 0.5*job->b13[i]*job->b21[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[0][6] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b21[i]*job->b24[i] - 1.0*job->b11[i]*job->b14[i] - 0.5*job->b21[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[0][7] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b11[i]*job->b24[i] + 0.5*NUMOD*job->b14[i]*job->b21[i] - 0.5*job->b14[i]*job->b21[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[1][0] += ( 0.5*EMOD*job->b11[i]*job->b21[i]*job->particles[i].v/(-NUMOD + 1) );
    job->elements[p].kku_element[1][1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*pow(job->b11[i], 2) - 0.5*pow(job->b11[i], 2) - 1.0*pow(job->b21[i], 2))/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[1][2] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b22[i] - 1.0*NUMOD*job->b12[i]*job->b21[i] - 0.5*job->b11[i]*job->b22[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[1][3] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b12[i] - 0.5*job->b11[i]*job->b12[i] - 1.0*job->b21[i]*job->b22[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[1][4] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b23[i] - 1.0*NUMOD*job->b13[i]*job->b21[i] - 0.5*job->b11[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[1][5] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b13[i] - 0.5*job->b11[i]*job->b13[i] - 1.0*job->b21[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[1][6] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b24[i] - 1.0*NUMOD*job->b14[i]*job->b21[i] - 0.5*job->b11[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[1][7] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b14[i] - 0.5*job->b11[i]*job->b14[i] - 1.0*job->b21[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[2][0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b21[i]*job->b22[i] - 1.0*job->b11[i]*job->b12[i] - 0.5*job->b21[i]*job->b22[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[2][1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b22[i] - 1.0*NUMOD*job->b12[i]*job->b21[i] - 0.5*job->b11[i]*job->b22[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[2][2] += ( EMOD*job->particles[i].v*(0.5*NUMOD*pow(job->b22[i], 2) - 1.0*pow(job->b12[i], 2) - 0.5*pow(job->b22[i], 2))/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[2][3] += ( 0.5*EMOD*job->b12[i]*job->b22[i]*job->particles[i].v/(-NUMOD + 1) );
    job->elements[p].kku_element[2][4] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b22[i]*job->b23[i] - 1.0*job->b12[i]*job->b13[i] - 0.5*job->b22[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[2][5] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b12[i]*job->b23[i] + 0.5*NUMOD*job->b13[i]*job->b22[i] - 0.5*job->b13[i]*job->b22[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[2][6] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b22[i]*job->b24[i] - 1.0*job->b12[i]*job->b14[i] - 0.5*job->b22[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[2][7] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b12[i]*job->b24[i] + 0.5*NUMOD*job->b14[i]*job->b22[i] - 0.5*job->b14[i]*job->b22[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[3][0] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b11[i]*job->b22[i] + 0.5*NUMOD*job->b12[i]*job->b21[i] - 0.5*job->b12[i]*job->b21[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[3][1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b12[i] - 0.5*job->b11[i]*job->b12[i] - 1.0*job->b21[i]*job->b22[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[3][2] += ( 0.5*EMOD*job->b12[i]*job->b22[i]*job->particles[i].v/(-NUMOD + 1) );
    job->elements[p].kku_element[3][3] += ( EMOD*job->particles[i].v*(0.5*NUMOD*pow(job->b12[i], 2) - 0.5*pow(job->b12[i], 2) - 1.0*pow(job->b22[i], 2))/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[3][4] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b12[i]*job->b23[i] - 1.0*NUMOD*job->b13[i]*job->b22[i] - 0.5*job->b12[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[3][5] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b12[i]*job->b13[i] - 0.5*job->b12[i]*job->b13[i] - 1.0*job->b22[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[3][6] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b12[i]*job->b24[i] - 1.0*NUMOD*job->b14[i]*job->b22[i] - 0.5*job->b12[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[3][7] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b12[i]*job->b14[i] - 0.5*job->b12[i]*job->b14[i] - 1.0*job->b22[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[4][0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b21[i]*job->b23[i] - 1.0*job->b11[i]*job->b13[i] - 0.5*job->b21[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[4][1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b23[i] - 1.0*NUMOD*job->b13[i]*job->b21[i] - 0.5*job->b11[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[4][2] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b22[i]*job->b23[i] - 1.0*job->b12[i]*job->b13[i] - 0.5*job->b22[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[4][3] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b12[i]*job->b23[i] - 1.0*NUMOD*job->b13[i]*job->b22[i] - 0.5*job->b12[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[4][4] += ( EMOD*job->particles[i].v*(0.5*NUMOD*pow(job->b23[i], 2) - 1.0*pow(job->b13[i], 2) - 0.5*pow(job->b23[i], 2))/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[4][5] += ( 0.5*EMOD*job->b13[i]*job->b23[i]*job->particles[i].v/(-NUMOD + 1) );
    job->elements[p].kku_element[4][6] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b23[i]*job->b24[i] - 1.0*job->b13[i]*job->b14[i] - 0.5*job->b23[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[4][7] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b13[i]*job->b24[i] + 0.5*NUMOD*job->b14[i]*job->b23[i] - 0.5*job->b14[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[5][0] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b11[i]*job->b23[i] + 0.5*NUMOD*job->b13[i]*job->b21[i] - 0.5*job->b13[i]*job->b21[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[5][1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b13[i] - 0.5*job->b11[i]*job->b13[i] - 1.0*job->b21[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[5][2] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b12[i]*job->b23[i] + 0.5*NUMOD*job->b13[i]*job->b22[i] - 0.5*job->b13[i]*job->b22[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[5][3] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b12[i]*job->b13[i] - 0.5*job->b12[i]*job->b13[i] - 1.0*job->b22[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[5][4] += ( 0.5*EMOD*job->b13[i]*job->b23[i]*job->particles[i].v/(-NUMOD + 1) );
    job->elements[p].kku_element[5][5] += ( EMOD*job->particles[i].v*(0.5*NUMOD*pow(job->b13[i], 2) - 0.5*pow(job->b13[i], 2) - 1.0*pow(job->b23[i], 2))/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[5][6] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b13[i]*job->b24[i] - 1.0*NUMOD*job->b14[i]*job->b23[i] - 0.5*job->b13[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[5][7] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b13[i]*job->b14[i] - 0.5*job->b13[i]*job->b14[i] - 1.0*job->b23[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[6][0] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b21[i]*job->b24[i] - 1.0*job->b11[i]*job->b14[i] - 0.5*job->b21[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[6][1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b24[i] - 1.0*NUMOD*job->b14[i]*job->b21[i] - 0.5*job->b11[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[6][2] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b22[i]*job->b24[i] - 1.0*job->b12[i]*job->b14[i] - 0.5*job->b22[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[6][3] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b12[i]*job->b24[i] - 1.0*NUMOD*job->b14[i]*job->b22[i] - 0.5*job->b12[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[6][4] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b23[i]*job->b24[i] - 1.0*job->b13[i]*job->b14[i] - 0.5*job->b23[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[6][5] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b13[i]*job->b24[i] - 1.0*NUMOD*job->b14[i]*job->b23[i] - 0.5*job->b13[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[6][6] += ( EMOD*job->particles[i].v*(0.5*NUMOD*pow(job->b24[i], 2) - 1.0*pow(job->b14[i], 2) - 0.5*pow(job->b24[i], 2))/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[6][7] += ( 0.5*EMOD*job->b14[i]*job->b24[i]*job->particles[i].v/(-NUMOD + 1) );
    job->elements[p].kku_element[7][0] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b11[i]*job->b24[i] + 0.5*NUMOD*job->b14[i]*job->b21[i] - 0.5*job->b14[i]*job->b21[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[7][1] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b11[i]*job->b14[i] - 0.5*job->b11[i]*job->b14[i] - 1.0*job->b21[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[7][2] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b12[i]*job->b24[i] + 0.5*NUMOD*job->b14[i]*job->b22[i] - 0.5*job->b14[i]*job->b22[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[7][3] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b12[i]*job->b14[i] - 0.5*job->b12[i]*job->b14[i] - 1.0*job->b22[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[7][4] += ( EMOD*job->particles[i].v*(-1.0*NUMOD*job->b13[i]*job->b24[i] + 0.5*NUMOD*job->b14[i]*job->b23[i] - 0.5*job->b14[i]*job->b23[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[7][5] += ( EMOD*job->particles[i].v*(0.5*NUMOD*job->b13[i]*job->b14[i] - 0.5*job->b13[i]*job->b14[i] - 1.0*job->b23[i]*job->b24[i])/(pow(NUMOD, 2) - 1) );
    job->elements[p].kku_element[7][6] += ( 0.5*EMOD*job->b14[i]*job->b24[i]*job->particles[i].v/(-NUMOD + 1) );
    job->elements[p].kku_element[7][7] += ( EMOD*job->particles[i].v*(0.5*NUMOD*pow(job->b14[i], 2) - 0.5*pow(job->b14[i], 2) - 1.0*pow(job->b24[i], 2))/(pow(NUMOD, 2) - 1) );

    /* --- K_geo --- Symmetric? True --- */
    job->elements[p].kku_element[0][0] += ( job->particles[i].v*(pow(job->b11[i], 2)*job->particles[i].sxx + 2*job->b11[i]*job->b21[i]*job->particles[i].sxy + pow(job->b21[i], 2)*job->particles[i].syy) );
    job->elements[p].kku_element[0][2] += ( job->particles[i].v*(job->b11[i]*job->b12[i]*job->particles[i].sxx + job->b11[i]*job->b22[i]*job->particles[i].sxy + job->b12[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b22[i]*job->particles[i].syy) );
    job->elements[p].kku_element[0][4] += ( job->particles[i].v*(job->b11[i]*job->b13[i]*job->particles[i].sxx + job->b11[i]*job->b23[i]*job->particles[i].sxy + job->b13[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b23[i]*job->particles[i].syy) );
    job->elements[p].kku_element[0][6] += ( job->particles[i].v*(job->b11[i]*job->b14[i]*job->particles[i].sxx + job->b11[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b24[i]*job->particles[i].syy) );
    job->elements[p].kku_element[1][1] += ( job->particles[i].v*(pow(job->b11[i], 2)*job->particles[i].sxx + 2*job->b11[i]*job->b21[i]*job->particles[i].sxy + pow(job->b21[i], 2)*job->particles[i].syy) );
    job->elements[p].kku_element[1][3] += ( job->particles[i].v*(job->b11[i]*job->b12[i]*job->particles[i].sxx + job->b11[i]*job->b22[i]*job->particles[i].sxy + job->b12[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b22[i]*job->particles[i].syy) );
    job->elements[p].kku_element[1][5] += ( job->particles[i].v*(job->b11[i]*job->b13[i]*job->particles[i].sxx + job->b11[i]*job->b23[i]*job->particles[i].sxy + job->b13[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b23[i]*job->particles[i].syy) );
    job->elements[p].kku_element[1][7] += ( job->particles[i].v*(job->b11[i]*job->b14[i]*job->particles[i].sxx + job->b11[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b24[i]*job->particles[i].syy) );
    job->elements[p].kku_element[2][0] += ( job->particles[i].v*(job->b11[i]*job->b12[i]*job->particles[i].sxx + job->b11[i]*job->b22[i]*job->particles[i].sxy + job->b12[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b22[i]*job->particles[i].syy) );
    job->elements[p].kku_element[2][2] += ( job->particles[i].v*(pow(job->b12[i], 2)*job->particles[i].sxx + 2*job->b12[i]*job->b22[i]*job->particles[i].sxy + pow(job->b22[i], 2)*job->particles[i].syy) );
    job->elements[p].kku_element[2][4] += ( job->particles[i].v*(job->b12[i]*job->b13[i]*job->particles[i].sxx + job->b12[i]*job->b23[i]*job->particles[i].sxy + job->b13[i]*job->b22[i]*job->particles[i].sxy + job->b22[i]*job->b23[i]*job->particles[i].syy) );
    job->elements[p].kku_element[2][6] += ( job->particles[i].v*(job->b12[i]*job->b14[i]*job->particles[i].sxx + job->b12[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b22[i]*job->particles[i].sxy + job->b22[i]*job->b24[i]*job->particles[i].syy) );
    job->elements[p].kku_element[3][1] += ( job->particles[i].v*(job->b11[i]*job->b12[i]*job->particles[i].sxx + job->b11[i]*job->b22[i]*job->particles[i].sxy + job->b12[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b22[i]*job->particles[i].syy) );
    job->elements[p].kku_element[3][3] += ( job->particles[i].v*(pow(job->b12[i], 2)*job->particles[i].sxx + 2*job->b12[i]*job->b22[i]*job->particles[i].sxy + pow(job->b22[i], 2)*job->particles[i].syy) );
    job->elements[p].kku_element[3][5] += ( job->particles[i].v*(job->b12[i]*job->b13[i]*job->particles[i].sxx + job->b12[i]*job->b23[i]*job->particles[i].sxy + job->b13[i]*job->b22[i]*job->particles[i].sxy + job->b22[i]*job->b23[i]*job->particles[i].syy) );
    job->elements[p].kku_element[3][7] += ( job->particles[i].v*(job->b12[i]*job->b14[i]*job->particles[i].sxx + job->b12[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b22[i]*job->particles[i].sxy + job->b22[i]*job->b24[i]*job->particles[i].syy) );
    job->elements[p].kku_element[4][0] += ( job->particles[i].v*(job->b11[i]*job->b13[i]*job->particles[i].sxx + job->b11[i]*job->b23[i]*job->particles[i].sxy + job->b13[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b23[i]*job->particles[i].syy) );
    job->elements[p].kku_element[4][2] += ( job->particles[i].v*(job->b12[i]*job->b13[i]*job->particles[i].sxx + job->b12[i]*job->b23[i]*job->particles[i].sxy + job->b13[i]*job->b22[i]*job->particles[i].sxy + job->b22[i]*job->b23[i]*job->particles[i].syy) );
    job->elements[p].kku_element[4][4] += ( job->particles[i].v*(pow(job->b13[i], 2)*job->particles[i].sxx + 2*job->b13[i]*job->b23[i]*job->particles[i].sxy + pow(job->b23[i], 2)*job->particles[i].syy) );
    job->elements[p].kku_element[4][6] += ( job->particles[i].v*(job->b13[i]*job->b14[i]*job->particles[i].sxx + job->b13[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b23[i]*job->particles[i].sxy + job->b23[i]*job->b24[i]*job->particles[i].syy) );
    job->elements[p].kku_element[5][1] += ( job->particles[i].v*(job->b11[i]*job->b13[i]*job->particles[i].sxx + job->b11[i]*job->b23[i]*job->particles[i].sxy + job->b13[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b23[i]*job->particles[i].syy) );
    job->elements[p].kku_element[5][3] += ( job->particles[i].v*(job->b12[i]*job->b13[i]*job->particles[i].sxx + job->b12[i]*job->b23[i]*job->particles[i].sxy + job->b13[i]*job->b22[i]*job->particles[i].sxy + job->b22[i]*job->b23[i]*job->particles[i].syy) );
    job->elements[p].kku_element[5][5] += ( job->particles[i].v*(pow(job->b13[i], 2)*job->particles[i].sxx + 2*job->b13[i]*job->b23[i]*job->particles[i].sxy + pow(job->b23[i], 2)*job->particles[i].syy) );
    job->elements[p].kku_element[5][7] += ( job->particles[i].v*(job->b13[i]*job->b14[i]*job->particles[i].sxx + job->b13[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b23[i]*job->particles[i].sxy + job->b23[i]*job->b24[i]*job->particles[i].syy) );
    job->elements[p].kku_element[6][0] += ( job->particles[i].v*(job->b11[i]*job->b14[i]*job->particles[i].sxx + job->b11[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b24[i]*job->particles[i].syy) );
    job->elements[p].kku_element[6][2] += ( job->particles[i].v*(job->b12[i]*job->b14[i]*job->particles[i].sxx + job->b12[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b22[i]*job->particles[i].sxy + job->b22[i]*job->b24[i]*job->particles[i].syy) );
    job->elements[p].kku_element[6][4] += ( job->particles[i].v*(job->b13[i]*job->b14[i]*job->particles[i].sxx + job->b13[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b23[i]*job->particles[i].sxy + job->b23[i]*job->b24[i]*job->particles[i].syy) );
    job->elements[p].kku_element[6][6] += ( job->particles[i].v*(pow(job->b14[i], 2)*job->particles[i].sxx + 2*job->b14[i]*job->b24[i]*job->particles[i].sxy + pow(job->b24[i], 2)*job->particles[i].syy) );
    job->elements[p].kku_element[7][1] += ( job->particles[i].v*(job->b11[i]*job->b14[i]*job->particles[i].sxx + job->b11[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b21[i]*job->particles[i].sxy + job->b21[i]*job->b24[i]*job->particles[i].syy) );
    job->elements[p].kku_element[7][3] += ( job->particles[i].v*(job->b12[i]*job->b14[i]*job->particles[i].sxx + job->b12[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b22[i]*job->particles[i].sxy + job->b22[i]*job->b24[i]*job->particles[i].syy) );
    job->elements[p].kku_element[7][5] += ( job->particles[i].v*(job->b13[i]*job->b14[i]*job->particles[i].sxx + job->b13[i]*job->b24[i]*job->particles[i].sxy + job->b14[i]*job->b23[i]*job->particles[i].sxy + job->b23[i]*job->b24[i]*job->particles[i].syy) );
    job->elements[p].kku_element[7][7] += ( job->particles[i].v*(pow(job->b14[i], 2)*job->particles[i].sxx + 2*job->b14[i]*job->b24[i]*job->particles[i].sxy + pow(job->b24[i], 2)*job->particles[i].syy) );


    #define DISPLACEMENT_DOFS 2

    #undef DISPLACEMENT_DOFS

    return;
}
/*----------------------------------------------------------------------------*/

