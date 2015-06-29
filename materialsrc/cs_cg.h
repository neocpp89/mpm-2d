/**
    \file cs_cg.h
    \author Sachith Dunatunga
    \date 06.05.2015

    A simple conjugate-gradient solver for CSparse.
*/
#include <suitesparse/cs.h>

double dot(const double * a, const double * b, size_t n);
int cs_cg(const cs *K, double *f, const double *u_0, double tol);
int cs_bicgstab(const cs *A, double *b, const double *x_0, double tol);
