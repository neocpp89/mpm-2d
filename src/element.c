/**
    \file element.c
    \author Sachith Dunatunga
    \date 04.06.12

    mpm_2d -- An implementation of the Material Point Method in 2D.
*/

/*---node_number_to_coords----------------------------------------------------*/
void node_number_to_coords(double *x, double *y, int num, int N, double h)
{
    const int i = num % N;
    const int j = num / N;

    *x = i*h;
    *y = j*h;

    return;
}
/*----------------------------------------------------------------------------*/

