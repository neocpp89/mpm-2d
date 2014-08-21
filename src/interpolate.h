/**
    \file interpolate.h
    \author Sachith Dunatunga
    \date 04.06.12

    Interpolation functions and gradient of interpolation functions.
*/
#ifndef __INTERPOLATE_H__
#define __INTERPOLATE_H__

void tent(double *h1, double *h2, double *h3, double *h4,
    double x_local, double y_local);

void grad_tent(double *b11, double *b12, double *b13, double *b14,
    double *b21, double *b22, double *b23, double *b24,
    double x_local, double y_local, double h);

void paraboloid(double *h1, double *h2, double *h3, double *h4,
    double *h5, double *h6, double *h7, double *h8, double *h9,
    double x_local, double y_local);

void grad_paraboloid(
    double *b11, double *b12, double *b13, double *b14, 
    double *b15, double *b16, double *b17, double *b18,
    double *b19,
    double *b21, double *b22, double *b23, double *b24,
    double *b25, double *b26, double *b27, double *b28,
    double *b29,
    double x_local, double y_local, double h);

inline double si(double xp, double xi, double h)
{
    double r = 0;
    double xl;
    
    xl = xp - xi;
    
    if (0 <= xl  && xl <= h) {
        r = 1.0 - (xl / h);
    } else if (-h <= xl && xl <= 0) {
        r = 1 + (xl / h);
    }

    return r;
}

inline double si2d(double xp, double yp, double xi, double yi, double h)
{
    return si(xp, xi, h)*si(yp, yi, h);
}

#endif

