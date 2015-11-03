/**
    \file interpolate.h
    \author Sachith Dunatunga
    \date 04.06.12

    Interpolation functions and gradient of interpolation functions.
*/
#ifndef __INTERPOLATE_H__
#define __INTERPOLATE_H__

/*
    Node order is

    4----3
    |    |
    |    |
    1----2

    (1 in lower left, increases in positive direction ccw).
*/

inline void tent(
    double * restrict h1, double * restrict h2, double * restrict h3, double * restrict h4,
    double x1, double y1,
    double x2, double y2,
    double x3, double y3,
    double x4, double y4,
    double hx, double hy,
    double xp, double yp
)
{
    const double da = (hx * hy);

    // *h1 = ((hx + x1 - xp) / hx) * ((hy + y1 - yp) / hy);
    *h1 = (hx + x1 - xp) * (hy + y1 - yp) / da;
    // *h2 = ((hx + xp - x2) / hx) * ((hy + y2 - yp) / hy);
    *h2 = (hx + xp - x2) * (hy + y2 - yp) / da;
    // *h3 = ((hx + xp - x3) / hx) * ((hy + yp - y3) / hy);
    *h3 = (hx + xp - x3) * (hy + yp - y3) / da;
    // *h4 = ((hx + x4 - xp) / hx) * ((hy + yp - y4) / hy);
    *h4 = (hx + x4 - xp) * (hy + yp - y4) / da;

    return;
}

inline void grad_tent(
    double * restrict b11, double * restrict b12, double * restrict b13, double * restrict b14,
    double * restrict b21, double * restrict b22, double * restrict b23, double * restrict b24,
    double x1, double y1,
    double x2, double y2,
    double x3, double y3,
    double x4, double y4,
    double hx, double hy,
    double xp, double yp
)
{
    const double da = (hx * hy);

    // *h1 = ((hx + x1 - xp) / hx) * ((hy + y1 - yp) / hy);
    // *h2 = ((hx + xp - x2) / hx) * ((hy + y2 - yp) / hy);
    // *h3 = ((hx + xp - x3) / hx) * ((hy + yp - y3) / hy);
    // *h4 = ((hx + x4 - xp) / hx) * ((hy + yp - y4) / hy);

    /*
        d/dx1
    */
    *b11 = - (hy + y1 - yp) / da;
    *b12 = (hy + y2 - yp) / da;
    *b13 = (hy + yp - y3) / da;
    *b14 = - (hy + yp - y4) / da;

    /*
        d/dx2
    */
    *b21 = - (hx + x1 - xp) / da;
    *b22 = - (hx + xp - x2) / da;
    *b23 = (hx + xp - x3) / da;
    *b24 = (hx + x4 - xp) / da;

    return;
}
#endif

