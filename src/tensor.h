/**
    \file tensor.h
    \author Sachith Dunatunga
    \date 04.06.12

    mpm_2d -- An implementation of the Material Point Method in 2D.
*/
#ifndef __TENSOR_H__
#define __TENSOR_H__

void zero_array(double *array, size_t numel);
void tensor_trace(double* restrict trA, const double* restrict A, size_t dim);
void tensor_trace3(double* restrict trA, const double* restrict A);
void tensor_add(double *C, const double *A, const double *B, size_t dim);
void tensor_add3(double *C, const double *A, const double *B);
void tensor_copy(double* restrict C, const double* restrict A, size_t dim);
void tensor_copy3(double* restrict C, const double* restrict A);
void tensor_decompose(double* restrict A_0, double* restrict c, const double* restrict A, size_t dim);
void tensor_decompose3(double* restrict A_0, double* restrict c, const double* restrict A);
void tensor_multiply(double* restrict C, const double* restrict A, const double* restrict B, size_t dim);
void tensor_multiply3(double* restrict C, const double* restrict A, const double* restrict B);
void tensor_multiply3_helper(double* restrict C, const double* restrict A, const double* restrict B);
void tensor_scale(double *A, double c, size_t dim);
void tensor_scale3(double *A, double c);
void tensor_sym(double * restrict C, const double * restrict A, size_t dim);
void tensor_sym3(double * restrict C, const double * restrict A);
void tensor_skw(double * restrict C, const double * restrict A, size_t dim);
void tensor_skw3(double * restrict C, const double * restrict A);
void tensor_contraction(double *c, double *A, double *B, size_t dim);
void tensor_contraction3(double *c, const double *A, const double *B);

#endif // __TENSOR_H__
