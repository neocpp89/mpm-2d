/**
    \file test_material.c
    \author Sachith Dunatunga
    \date 26.06.2014

    Test a material model with given velocity gradients.
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

void calculate_stress(job_t *job);

double mu(particle_t *p)
{
    double pr = -0.5 * (p->sxx + p->syy);
    double s0xx = p->sxx + pr;
    double s0xy = p->sxy;
    double s0yy = p->syy + pr;
    double tau = sqrt(0.5*(s0xx*s0xx + 2*s0xy*s0xy + s0yy*s0yy));

    if (pr <= 0) {
        return 0;
    }
    return tau/pr;
}

void set_velocity_gradient(particle_t *p, double t)
{
    double sh = 0.5;
    double L[4] = { 0, 0, 0, 0 };
    int i;

    if (t >= 0 && t < 0.25) {
        L[0] = 0;
/*        L[1] = 0;*/
        L[1] = sh;
        L[2] = 0;
        L[3] = 0;
/*        L[3] = -sh;*/
    } else if (t >= 0.25 && t < 0.50) {
        L[0] = -1*(0.125 - fabs(t - 0.375))/0.125;
        L[1] = sh;
        L[2] = 0;
        L[3] = -1*(0.125 - fabs(t - 0.375))/0.125;
    } else if (t >= 0.50 && t < 0.75) {
        L[0] = 0;
        L[1] = sh;
        L[2] = 0;
        L[3] = 0;
    } else if (t >= 0.75 && t < 1) {
        L[0] = 1*(0.125 - fabs(t - 0.875))/0.125;
        L[1] = sh;
        L[2] = 0;
        L[3] = 1*(0.125 - fabs(t - 0.875))/0.125;
    } else {
        L[0] = 0;
        L[1] = sh;
        L[2] = 0;
        L[3] = 0;
    }

    for (i = 0; i < 4; i++) {
        L[i] *= 1e-1;
    }

    p->exx_t = L[0];
    p->exy_t = 0.5 * (L[1] + L[2]);
    p->wxy_t = 0.5 * (L[1] - L[2]);
    p->eyy_t = L[3];

    return;
}

void update_density(particle_t *p, double dt)
{
    p->v *= exp (dt * (p->exx_t + p->eyy_t));
    return;
}

int main(int argc, char **argv)
{
    job_t testjob;
/*    double testprops[2] = { 72.0e3/7.0, 2.0/7.0 };*/
    double testprops[2] = { 1e4, 0.3 };
    int a[1] = {1};
    double t_stop = 1.25;
    particle_t p;

    p.sxx = -1;
    p.sxy = 0;
    p.syy = -1;
    p.m = 1500;
    p.v = 1;

    testjob.num_particles = 1;
    testjob.active = a;
    testjob.particles = &p;
    testjob.material.num_fp64_props = 2;
    testjob.material.fp64_props = testprops;
    testjob.t = 0;
    testjob.dt = 1e-3;


    if (argc > 1) {
        testjob.dt = atof(argv[1]);
    }

    if (testjob.dt <= 0) {
        testjob.dt = 1e-3;
    }

    FILE *fp = fopen("mattest_output.csv", "w");

    material_init(&testjob);
    fprintf(fp, "%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
        testjob.t,
        p.exx_t, p.exy_t + p.wxy_t, p.exy_t - p.wxy_t, p.eyy_t,
        p.sxx, p.sxy, p.syy, mu(&p),
        p.state[10] /* gammadotp */
    );

    while(testjob.t < t_stop) {
        set_velocity_gradient(&p, testjob.t);
        calculate_stress(&testjob);
        update_density(&p, testjob.dt);
        testjob.t += testjob.dt;
        fprintf(fp, "%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
            testjob.t,
            p.exx_t, p.exy_t + p.wxy_t, p.exy_t - p.wxy_t, p.eyy_t,
            p.sxx, p.sxy, p.syy, mu(&p),
            p.state[10] /* gammadotp */
        );
    }

    fclose(fp);

    return 0;
}

