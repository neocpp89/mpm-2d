/**
    \file interpolate.c
    \author Sachith Dunatunga
    \date 04.06.12

    mpm_2d -- An implementation of the Material Point Method in 2D.
*/
#include <stdio.h>
#include "interpolate.h"

#define SF_WARNING4(tok) \
    if (*tok > 1 || *tok < 0) { \
        fprintf(stderr, "warning, " #tok " = %g\n", *tok); \
    }

#define SF_WARNING(tok) \
    if (*tok > 1 || *tok < -1) { \
        fprintf(stderr, "warning, " #tok " = %g\n", *tok); \
    }

/*---tent---------------------------------------------------------------------*/
void tent(double *h1, double *h2, double *h3, double *h4,
    double x_local, double y_local)
{
    *h1 = (1 - x_local) * (1 - y_local);
    *h2 = (x_local) * (1 - y_local);
    *h3 = (x_local) * (y_local);
    *h4 = (1 - x_local) * (y_local);

    SF_WARNING4(h1);
    SF_WARNING4(h2);
    SF_WARNING4(h3);
    SF_WARNING4(h4);

    return;
}
/*----------------------------------------------------------------------------*/

/*---paraboloid---------------------------------------------------------------*/
void paraboloid(double *h1, double *h2, double *h3, double *h4,
    double *h5, double *h6, double *h7, double *h8, double *h9,
    double x_local, double y_local)
{
    /*
        1,2,3,4 are corners from bottom-left proceeding in positive angular
        direction (counter-clockwise).

        5,6,7,8 are edge nodes from bottom edge proceeding in positive angular
        direction (counter-clockwise).

        9 is the center node. The origin of x_local and y_local is at node 9 in
        a square element, and x_local and y_local are in [-1,1].
    */

    *h1 = 0.25 * (1 - x_local) * (1 - y_local) * x_local * y_local;
    *h2 = -0.25 * (1 + x_local) * (1 - y_local) * x_local * y_local;
    *h3 = 0.25 * (1 + x_local) * (1 + y_local) * x_local * y_local;
    *h4 = -0.25 * (1 - x_local) * (1 + y_local) * x_local * y_local;

    *h5 = -0.5 * (1 - x_local * x_local) * (1 - y_local) * y_local;
    *h6 = 0.5 * (1 - y_local * y_local) * (1 + x_local) * x_local;
    *h7 = 0.5 * (1 - x_local * x_local) * (1 + y_local) * y_local;
    *h8 = -0.5 * (1 - y_local * y_local) * (1 - x_local) * x_local;

    *h9 = (1 - x_local * x_local) * (1 - y_local * y_local);

    SF_WARNING(h1);
    SF_WARNING(h2);
    SF_WARNING(h3);
    SF_WARNING(h4);
    SF_WARNING(h5);
    SF_WARNING(h6);
    SF_WARNING(h7);
    SF_WARNING(h8);
    SF_WARNING(h9);

    return;
}
/*----------------------------------------------------------------------------*/

/*---grad_tent----------------------------------------------------------------*/
void grad_tent(double *b11, double *b12, double *b13, double *b14,
    double *b21, double *b22, double *b23, double *b24,
    double x_local, double y_local, double h)
{
    double dxl_dx = (1.0f / h);
    double dyl_dy = (1.0f / h);

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

/*---grad_paraboloid----------------------------------------------------------*/
void grad_paraboloid(
    double *b11, double *b12, double *b13, double *b14, 
    double *b15, double *b16, double *b17, double *b18,
    double *b19,
    double *b21, double *b22, double *b23, double *b24,
    double *b25, double *b26, double *b27, double *b28,
    double *b29,
    double x_local, double y_local, double h)
{
    /*
        *h1 = 0.25 * (1 - x_local) * (1 - y_local) * x_local * y_local;
        *h2 = -0.25 * (1 + x_local) * (1 - y_local) * x_local * y_local;
        *h3 = 0.25 * (1 + x_local) * (1 + y_local) * x_local * y_local;
        *h4 = -0.25 * (1 - x_local) * (1 + y_local) * x_local * y_local;

        *h5 = -0.5 * (1 - x_local * x_local) * (1 - y_local) * y_local;
        *h6 = 0.5 * (1 - y_local * y_local) * (1 + x_local) * x_local;
        *h7 = 0.5 * (1 - x_local * x_local) * (1 + y_local) * y_local;
        *h8 = -0.5 * (1 - y_local * y_local) * (1 - x_local) * x_local;

        *h9 = (1 - x_local * x_local) * (1 - y_local * y_local);

        h is spacing between grid points NOT size of element!
        size of element is then given by 2h x 2h.
    */
    double dxl_dx = (1.0 / h);
    double dyl_dy = (1.0 / h);

    /*
        d/dx1
    */
    *b11 = 0.25 * (1 - y_local) * y_local * (1 - 2 * x_local) * dxl_dx;
    *b12 = -0.25 * (1 - y_local) * y_local * (1 + 2 * x_local) * dxl_dx;
    *b13 = 0.25 * (1 + y_local) * y_local * (1 + 2 * x_local) * dxl_dx;
    *b14 = -0.25 * (1 + y_local) * y_local * (1 - 2 * x_local) *dxl_dx;

    *b15 = x_local * (1 - y_local) * y_local * dxl_dx;
    *b16 = 0.5 * (1 - y_local * y_local) * (1 + 2 * x_local) * dxl_dx;
    *b17 = -x_local * (1 + y_local) * y_local * dxl_dx;
    *b18 = -0.5 * (1 - y_local * y_local) * (1 - 2 * x_local) * dxl_dx;

    *b19 = -2 * x_local * (1 - y_local * y_local) * dxl_dx;

    /*
        d/dx2
    */
    *b21 = 0.25 * (1 - x_local) * x_local * (1 - 2 * y_local) * dyl_dy;
    *b22 = -0.25 * (1 + x_local) * x_local * (1 - 2 * y_local) * dyl_dy;
    *b23 = 0.25 * (1 + x_local) * x_local * (1 + 2 * y_local) * dyl_dy;
    *b24 = -0.25 * (1 - x_local) * x_local * (1 + 2 * y_local) * dyl_dy;
    
    *b25 = -0.5 * (1 - x_local * x_local) * (1 - 2 * y_local) * dyl_dy;
    *b26 = -y_local * (1 + x_local) * x_local * dyl_dy;
    *b27 = 0.5 * (1 - x_local * x_local) * (1 + 2 * y_local) * dyl_dy;
    *b28 = y_local * (1 - x_local) * x_local * dyl_dy;

    *b29 = -2 * y_local * (1 - x_local * x_local) * dyl_dy;

    return;
}
/*----------------------------------------------------------------------------*/

