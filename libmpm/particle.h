/**
    \file particle.h
    \author Sachith Dunatunga
    \date 04.06.12

    Contains the structure for particles in MPM.
*/
#include <stddef.h>

#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#define DEPVAR 11

/* shape function related indicies in arrays */
#define S_XIDX 0
#define S_YIDX 1
#define S_ZIDX 2

/* dimension of the stress/strain tensors */
#define NDIM 3

/* indicies for stress/strain tensors */
#define XX (NDIM * S_XIDX + S_XIDX)
#define XY (NDIM * S_XIDX + S_YIDX)
#define XZ (NDIM * S_XIDX + S_ZIDX)
#define YX (NDIM * S_YIDX + S_XIDX)
#define YY (NDIM * S_YIDX + S_YIDX)
#define YZ (NDIM * S_YIDX + S_ZIDX)
#define ZX (NDIM * S_ZIDX + S_XIDX)
#define ZY (NDIM * S_ZIDX + S_YIDX)
#define ZZ (NDIM * S_ZIDX + S_ZIDX)

typedef struct particle_s {
    /* Position */
    double x;
    double y;

    /* Volume */
    double v;

    /* Initial volume */
    double v0;

    /* Shapefunctions (follows same nodal numbering as element) */
    double s[4];

    /* Mass */
    double m;

    /* Velocity */
    double x_t;
    double y_t;

    /* Acceleration */
    double x_tt;
    double y_tt;

    /* Body forces */
    double bx;
    double by;

    /* Stress */
    double sxx;
    double sxy;
    double syy;

    /* full 3D stress tensor */
    double T[NDIM*NDIM];

    /* Strain rate */
    double exx_t;
    double exy_t;
    double eyy_t;
    double wxy_t;

    /* full 3D velocity gradient tensor */
    double L[NDIM*NDIM];

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
//    int active;

    /* Material Type */
    int material;

    /* Gradient of Shapefunctions */
    double grad_s[4][2];

    /* which element is the particle in? */
    int element;

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
    
    /* unique id, so particles can be tracked between frames */
    size_t id;
} particle_t;

/*
    Convert global coordinates to local coordinates of a square element with
    size h by h and bottom left corner at (x_ref, y_ref).
*/
void global_to_local_coords(double *x_local, double *y_local, 
    double x, double y, 
    double x_ref, double y_ref, double hx, double hy);
#endif

