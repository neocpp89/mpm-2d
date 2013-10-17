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
void calculate_stress_linear_elastic_novisc(job_t *job)
{
    double txx;
    double txy;
    double tyy;

    double dsjxx;
    double dsjxy;
    double dsjyy;

    int i;

    for (i = 0; i < job->num_particles; i++) {
        dsjxx = job->dt * (EMOD / (1 - NUMOD*NUMOD)) * ((job->particles[i].exx_t) + NUMOD * (job->particles[i].eyy_t));
        dsjxy = job->dt * (EMOD / (2 *(1 + NUMOD))) * (job->particles[i].exy_t);
        dsjyy = job->dt * (EMOD / (1 - NUMOD*NUMOD)) * ((job->particles[i].eyy_t) + NUMOD * (job->particles[i].exx_t));

        dsjxx -= 2 * job->dt * job->particles[i].wxy_t * job->particles[i].sxy;
        dsjxy += job->dt * job->particles[i].wxy_t * (job->particles[i].sxx - job->particles[i].syy);
        dsjyy += 2 * job->dt * job->particles[i].wxy_t * job->particles[i].sxy;

/*        dsjxx += 1e0 * (jp(exx_t) + jp(eyy_t));*/
/*        dsjyy += 1e0 * (jp(eyy_t) + jp(exx_t));*/

        txx = job->particles[i].sxx + dsjxx;
        txy = job->particles[i].sxy + dsjxy;
        tyy = job->particles[i].syy + dsjyy;

        job->particles[i].sxx = txx;
        job->particles[i].sxy = txy;
        job->particles[i].syy = tyy;
    }

    return;
}
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
void calculate_stress_dp_shearonly(job_t *job)
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

    double const c = 0;
    int i;

    for (i = 0; i < job->num_particles; i++) {
        dsjxx = job->dt * (EMOD / (1 - NUMOD*NUMOD)) * ((job->particles[i].exx_t) + NUMOD * (job->particles[i].eyy_t));
        dsjxy = job->dt * (EMOD / (2 *(1 + NUMOD))) * (job->particles[i].exy_t);
        dsjyy = job->dt * (EMOD / (1 - NUMOD*NUMOD)) * ((job->particles[i].eyy_t) + NUMOD * (job->particles[i].exx_t));

        dsjxx -= 2 * job->dt * job->particles[i].wxy_t * job->particles[i].sxy;
        dsjxy += job->dt * job->particles[i].wxy_t * (job->particles[i].sxx - job->particles[i].syy);
        dsjyy += 2 * job->dt * job->particles[i].wxy_t * job->particles[i].sxy;

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
            dp_xx = 1e-6 * f * t0xx;
            dp_xy = 1e-6 * f * t0xy;
            dp_yy = 1e-6 * f * t0yy;
            job->particles[i].sxx = txx - job->dt * (EMOD / (1 - NUMOD*NUMOD)) * ((dp_xx) + NUMOD * (dp_yy));
            job->particles[i].sxy = txy - job->dt * (EMOD / (2 *(1 + NUMOD))) * (dp_xy);
            job->particles[i].syy = tyy - job->dt * (EMOD / (1 - NUMOD*NUMOD)) * ((dp_yy) + NUMOD * (dp_xx));
            gammap = sqrt(Epxx*Epxx + 2*Epxy*Epxy + Epyy*Epyy);
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void calculate_stress_dp_shearonly_indep(job_t *job)
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
    double q_adj;

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

/*        dsjxx += 2 * job->dt * job->particles[i].wxy_t * job->particles[i].sxy;*/
/*        dsjxy -= job->dt * job->particles[i].wxy_t * (job->particles[i].sxx - job->particles[i].syy);*/
/*        dsjyy -= 2 * job->dt * job->particles[i].wxy_t * job->particles[i].sxy;*/

/*        dsjxx += 1e2 * (jp(exx_t) + jp(eyy_t));*/
/*        dsjyy += 1e2 * (jp(eyy_t) + jp(exx_t));*/

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
/*            printf("particle %d is rigid.\n", i);*/
            continue;
        }

        if ((job->particles[i].m / job->particles[i].v) < 1485.0f) {
            density_flag = 1;
/*            printf("%4d: density %lf\n", i, (job->particles[i].m / job->particles[i].v));*/
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
            job->particles[i].color = 1;
        } else if (f >= 0 && pm > -c/m) {
            Epxx += (f / qm) * t0xx / (2 * G);
            Epxy += (f / qm) * t0xy / (2 * G);
            Epyy += (f / qm) * t0yy / (2 * G);
            q_adj = m*pm + c;
            job->particles[i].sxx = (q_adj / qm) * t0xx - pm;
            job->particles[i].sxy = (q_adj / qm) * t0xy;
            job->particles[i].syy = (q_adj / qm) * t0yy - pm;
            job->particles[i].color = 2;
            gammap = sqrt(Epxx*Epxx + 2*Epxy*Epxy + Epyy*Epyy);
        } else if (pm <= -c/m) {
            job->particles[i].sxx = -0.5 * c / m;
            job->particles[i].sxy = 0;
            job->particles[i].syy = -0.5 * c / m;
            job->particles[i].color = 3;
        } else {
            fprintf(stderr, "u");
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void calculate_stress_vm_indep(job_t *job)
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
    double q_adj;

    double dsjxx;
    double dsjxy;
    double dsjyy;

    int i;

    double q0 = 1e4;

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

        txx = job->particles[i].sxx + dsjxx;
        txy = job->particles[i].sxy + dsjxy;
        tyy = job->particles[i].syy + dsjyy;
        pm = -0.5f*(txx + tyy);
        t0xx = txx + pm;
        t0xy = txy;
        t0yy = tyy + pm;
        qm = sqrt(0.5f*(t0xx*t0xx + 2*t0xy*t0xy + t0yy*t0yy));

        f = qm - q0;

        if(job->particles[i].material == M_RIGID) {
            job->particles[i].sxx = txx;
            job->particles[i].sxy = txy;
            job->particles[i].syy = tyy;
/*            printf("particle %d is rigid.\n", i);*/
            continue;
        }

        if (f < 0) {
            job->particles[i].sxx = txx;
            job->particles[i].sxy = txy;
            job->particles[i].syy = tyy;
            job->particles[i].color = 1;
        } else if (f >= 0) {
            Epxx += (f / qm) * t0xx / (2 * G);
            Epxy += (f / qm) * t0xy / (2 * G);
            Epyy += (f / qm) * t0yy / (2 * G);
            q_adj = q0;
            job->particles[i].sxx = (q_adj / qm) * t0xx - pm;
            job->particles[i].sxy = (q_adj / qm) * t0xy;
            job->particles[i].syy = (q_adj / qm) * t0yy - pm;
            job->particles[i].color = 2;
            gammap = sqrt(Epxx*Epxx + 2*Epxy*Epxy + Epyy*Epyy);
        } else {
            fprintf(stderr, "undetermined case\n");
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

        job->particles[i].sxx = job->particles[i].sxx + dsjxx;
        job->particles[i].sxy = job->particles[i].sxy + dsjxy;
        job->particles[i].syy = job->particles[i].syy + dsjyy;
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void calculate_stress_dp_dilate(job_t *job)
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
    double a;
    double q_adj;
    double p_adj;

    double Ef_xx;
    double Ef_xy;
    double Ef_yy;

    double dsjxx;
    double dsjxy;
    double dsjyy;

    double density_flag;

    double const c = 1e-2;
    int i;

    for (i = 0; i < job->num_particles; i++) {
        dsjxx = job->dt * (EMOD / (1 - NUMOD*NUMOD)) * ((job->particles[i].exx_t) + NUMOD * (job->particles[i].eyy_t));
        dsjxy = job->dt * (EMOD / (2 *(1 + NUMOD))) * (job->particles[i].exy_t);
        dsjyy = job->dt * (EMOD / (1 - NUMOD*NUMOD)) * ((job->particles[i].eyy_t) + NUMOD * (job->particles[i].exx_t));

        dsjxx -= 2 * job->dt * job->particles[i].wxy_t * job->particles[i].sxy;
        dsjxy += job->dt * job->particles[i].wxy_t * (job->particles[i].sxx - job->particles[i].syy);
        dsjyy += 2 * job->dt * job->particles[i].wxy_t * job->particles[i].sxy;

/*        job->particles[i].sxx += 1e3 * (jp(exx_t) + jp(eyy_t));*/
/*        job->particles[i].syy += 1e3 * (jp(eyy_t) + jp(exx_t));*/
/*        job->particles[i].sxx += SVISC * jp(exx_t);*/
/*        job->particles[i].sxy += SVISC * jp(exy_t);*/
/*        job->particles[i].syy += SVISC * jp(eyy_t);*/

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

        if ((job->particles[i].m / job->particles[i].v) < 1485.0f) {
            density_flag = 1;
/*            printf("%4d: density %lf\n", i, (job->particles[i].m / job->particles[i].v));*/
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
            job->particles[i].color = 1;
        } else if (f >= 0 && pm > -c/m) {
            Epxx += (f / qm) * t0xx / (2 * G);
            Epxy += (f / qm) * t0xy / (2 * G);
            Epyy += (f / qm) * t0yy / (2 * G);
            a = qm + pm / m;
            q_adj = (a * m * m + c) / (m * m + 1.0f);
            p_adj = m * (a - c) / (m * m + 1.0f);
            job->particles[i].sxx = (q_adj / qm) * t0xx - p_adj;
            job->particles[i].sxy = (q_adj / qm) * t0xy;
            job->particles[i].syy = (q_adj / qm) * t0yy - p_adj;
/*            q_adj = m*pm + c;*/
/*            job->particles[i].sxx = (q_adj / qm) * t0xx - pm;*/
/*            job->particles[i].sxy = (q_adj / qm) * t0xy;*/
/*            job->particles[i].syy = (q_adj / qm) * t0yy - pm;*/
            job->particles[i].color = 2;
            gammap = sqrt(Epxx*Epxx + 2*Epxy*Epxy + Epyy*Epyy);
        } else if (pm <= -c/m) {
            job->particles[i].sxx = -0.5 * c / m;
            job->particles[i].sxy = 0;
            job->particles[i].syy = -0.5 * c / m;
            job->particles[i].color = 3;
        } else {
            fprintf(stderr, "undetermined case\n");
        }

        Ef_xx = 0.5f * (jp(Fxx)*jp(Fxx) + jp(Fxy)*jp(Fxy) - 1);
        Ef_xy = 0.5f * (jp(Fxx)*jp(Fyx) + jp(Fxy)*jp(Fyy));
        Ef_yy = 0.5f * (jp(Fyx)*jp(Fyx) + jp(Fyy)*jp(Fyy) - 1);

        Ef_mag =  sqrt(Ef_xx*Ef_xx + 2*Ef_xy*Ef_xy + Ef_yy*Ef_yy);

/*        job->particles[i].sxx += 1e5 * jp(exx_t);*/
/*        job->particles[i].sxy += 1e5 * jp(exy_t);*/
/*        job->particles[i].syy += 1e5 * jp(eyy_t);*/

/*        job->particles[i].color = gammap;*/
/*        job->particles[i].color = job->particles[i].syy;*/
    }

    return;
}
/*----------------------------------------------------------------------------*/

