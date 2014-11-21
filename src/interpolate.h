/**
    \file interpolate.h
    \author Sachith Dunatunga
    \date 04.06.12

    Interpolation functions and gradient of interpolation functions.
*/
#ifndef __INTERPOLATE_H__
#define __INTERPOLATE_H__

void tent(double * restrict h1, double * restrict h2, double * restrict h3, double * restrict h4,
    double x_local, double y_local);

void grad_tent(double * restrict b11, double * restrict b12, double * restrict b13, double * restrict b14,
    double * restrict b21, double * restrict b22, double * restrict b23, double * restrict b24,
    double x_local, double y_local, double h);

#endif

