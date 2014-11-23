/**
    \file interpolate.c
    \author Sachith Dunatunga
    \date 04.06.12

    mpm_2d -- An implementation of the Material Point Method in 2D.
*/
#include <stdio.h>
#include "interpolate.h"

#define SF_WARNING4(tok,xl,yl) \
    if (*tok > 1 || *tok < 0) { \
        fprintf(stderr, "warning, " #tok " = %g", *tok); \
        fprintf(stderr, ", (%g, %g)\n", xl, yl); \
    }

#define SF_WARNING(tok) \
    if (*tok > 1 || *tok < -1) { \
        fprintf(stderr, "warning, " #tok " = %g\n", *tok); \
    }

/*---tent---------------------------------------------------------------------*/
void tent(double * restrict h1, double * restrict h2, double * restrict h3, double * restrict h4,
    double x_local, double y_local)
{
    *h1 = (1 - x_local) * (1 - y_local);
    *h2 = (x_local) * (1 - y_local);
    *h3 = (x_local) * (y_local);
    *h4 = (1 - x_local) * (y_local);

    SF_WARNING4(h1, x_local, y_local);
    SF_WARNING4(h2, x_local, y_local);
    SF_WARNING4(h3, x_local, y_local);
    SF_WARNING4(h4, x_local, y_local);

    return;
}
/*----------------------------------------------------------------------------*/

/*---grad_tent----------------------------------------------------------------*/
void grad_tent(double * restrict b11, double * restrict b12, double * restrict b13, double * restrict b14,
    double * restrict b21, double * restrict b22, double * restrict b23, double * restrict b24,
    double x_local, double y_local, double h)
{
    double dxl_dx = (1.0 / h);
    double dyl_dy = (1.0 / h);

    /*
        d/dx1
    */
    *b11 = -(1 - y_local) * dxl_dx;
    *b12 = (1 - y_local) * dxl_dx;
    *b13 = (y_local) * dxl_dx;
    *b14 = -(y_local) * dxl_dx;

    /*
        d/dx2
    */
    *b21 = -(1 - x_local) * dyl_dy;
    *b22 = -(x_local) * dyl_dy;
    *b23 = (x_local) * dyl_dy;
    *b24 = (1 - x_local) * dyl_dy;

    return;
}
/*----------------------------------------------------------------------------*/

