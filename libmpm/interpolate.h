/**
    \file interpolate.h
    \author Sachith Dunatunga
    \date 04.06.12

    Interpolation functions and gradient of interpolation functions.
*/
#ifndef __INTERPOLATE_H__
#define __INTERPOLATE_H__

#if 0
void tent(double * restrict h1, double * restrict h2, double * restrict h3, double * restrict h4,
    double x_local, double y_local);

void grad_tent(double * restrict b11, double * restrict b12, double * restrict b13, double * restrict b14,
    double * restrict b21, double * restrict b22, double * restrict b23, double * restrict b24,
    double x_local, double y_local, double hx, double hy, double Lx, double Ly);
#endif

void tent(
    double * restrict h1, double * restrict h2, double * restrict h3, double * restrict h4,
    double x1, double y1,
    double x2, double y2,
    double x3, double y3,
    double x4, double y4,
    double hx, double hy,
    double xp, double yp
);

void grad_tent(
    double * restrict b11, double * restrict b12, double * restrict b13, double * restrict b14,
    double * restrict b21, double * restrict b22, double * restrict b23, double * restrict b24,
    double x1, double y1,
    double x2, double y2,
    double x3, double y3,
    double x4, double y4,
    double hx, double hy,
    double xp, double yp
);

#endif

