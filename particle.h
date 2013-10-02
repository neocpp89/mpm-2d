/**
    \file particle.h
    \author Sachith Dunatunga
    \date 04.06.12

    Contains the structure for particles in MPM.
*/
#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#define DEPVAR 11

typedef struct particle_s {
    /* Position */
    double x;
    double y;

    /* Local Coordinates */
    double xl;
    double yl;

    /* Velocity */
    double x_t;
    double y_t;

    /* Mass */
    double m;

    /* Volume */
    double v;

    /* Initial volume */
    double v0;

    /* Stress */
    double sxx;
    double sxy;
    double syy;

    /* Strain rate */
    double exx_t;
    double exy_t;
    double eyy_t;
    double wxy_t;

    /* Body forces */
    double bx;
    double by;

    /* Deformation gradient tensor */
    double Fxx;
    double Fxy;
    double Fyx;
    double Fyy;

    /* Displacements */
    double ux;
    double uy;

    /* Color used by splot visualization */
    double color;

    /* State Variables (for constitutive law) */
    double state[DEPVAR];

    /* Flag particle as active or not (used for discharge problems). */
    int active;

    /* Material Type */
    int material;

#if 0
    /* Shapefunctions */
    double h[9];

    /* Gradient of Shapefunctions */
    double b1[9];
    double b2[9];
#endif

    /* acceleration */
    double x_tt;
    double y_tt;

    /* held stress (for each implicit step) */
    double real_sxx;
    double real_sxy;
    double real_syy;
    double real_state[DEPVAR];

    /* inital 'vertical' and 'horizontal' vectors.
        initial volume should be given by magnitude (r_v cross r_h). */
    double r1_initial[2];
    double r2_initial[2];

    /* current vertical and horizontal vectors. */
    double r1[2];
    double r2[2];

    /* matrix of corner positions c[corner#][x or y] */
    double corners[4][2];

    /* which elements are the corners in? */
    int corner_elements[4];

} particle_t;

/*
    Convert global coordinates to local coordinates of a square element with
    size h by h and bottom left corner at (x_ref, y_ref).
*/
void global_to_local_coords(double *x_local, double *y_local, 
    double x, double y, 
    double x_ref, double y_ref, double h);
#endif

