/**
    \file element.c
    \author Sachith Dunatunga
    \date 04.06.12

    mpm_2d -- An implementation of the Material Point Method in 2D.
*/
#include "element.h"
#include "particle.h"

/*---node_number_to_coords----------------------------------------------------*/
void inline node_number_to_coords(double *x, double *y, int num, int N, double h)
{
    int i;
    int j;

    i = num % N;
    j = num / N;

    *x = i*h;
    *y = j*h;

    return;
}
/*----------------------------------------------------------------------------*/

