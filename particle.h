/**
    \file particle.h
    \author Sachith Dunatunga
    \date 04.06.12

    Contains the structure for particles in MPM.
*/
#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#define DEPVAR 11

/* shape function related indicies in arrays */
#define S_XIDX 0
#define S_YIDX 1

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

    /* Shapefunctions (follows same nodal numbering as element) */
    double s[4];

    /* Gradient of Shapefunctions */
    double grad_s[4][2];

    /* which element is the particle in? */
    int element;

    /* acceleration */
    double x_tt;
    double y_tt;

    /* held stress (for each implicit step) */
    double real_sxx;
    double real_sxy;
    double real_syy;
    double real_state[DEPVAR];

    /* inital 'horizontal' and 'vertical' vectors.
        initial volume should be given by magnitude (r_v cross r_h). */
    double r1_initial[2];
    double r2_initial[2];

    /* current vertical and horizontal vectors. */
    double r1[2];
    double r2[2];

    /* matrix of corner positions c[corner#][x or y] */
    double corners[4][2];

    /* matrix of corner positions in local coordinates in appropriate element */
    double cornersl[4][2];

    /* shapefunctions for corner points sc[corner#][node#] */
    double sc[4][4];

    /* gradient of shapefunctions for corner points sc[corner#][node#][x or y] */
    double grad_sc[4][4][2];

    /* which elements are the corners in? */
    int corner_elements[4];

    /* material specific data structure, can be defined in material_init */
    void *material_data;
} particle_t;

/*
    Convert global coordinates to local coordinates of a square element with
    size h by h and bottom left corner at (x_ref, y_ref).
*/
void global_to_local_coords(double *x_local, double *y_local, 
    double x, double y, 
    double x_ref, double y_ref, double h);
#endif

