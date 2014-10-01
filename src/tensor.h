/**
    \file tensor.h
    \author Sachith Dunatunga
    \date 04.06.12

    mpm_2d -- An implementation of the Material Point Method in 2D.
*/
#ifndef __TENSOR_H__
#define __TENSOR_H__

void zero_array(double *array, size_t numel);
void tensor_trace(double *trA, const double *A, size_t dim);
void tensor_add(double *C, const double *A, const double *B, size_t dim);
void tensor_copy(double *C, const double *A, size_t dim);
void tensor_decompose_deviatoric_and_spherical(double *A_0, double *c, double *A, size_t dim);

void tensor_add3(double* restrict C, const double * restrict A, const double * restrict B);
void tensor_copy3(double* restrict C, const double * restrict A);
void tensor_trace3(double * restrict trA, const double * restrict A);

void tensor_multiply(double* restrict C, const double* restrict A, const double* restrict B, size_t dim);
void tensor_multiply3(double* restrict C, const double* restrict A, const double* restrict B);
void tensor_multiply3_helper(double* restrict C, const double* restrict A, const double* restrict B);

#endif // __TENSOR_H__
