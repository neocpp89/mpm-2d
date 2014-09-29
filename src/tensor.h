/**
    \file tensor.h
    \author Sachith Dunatunga
    \date 04.06.12

    mpm_2d -- An implementation of the Material Point Method in 2D.
*/
#ifndef __TENSOR_H__
#define __TENSOR_H__

void zero_array(double *array, size_t numel);
void tensor_trace(double *trA, double *A, size_t dim);
void tensor_add(double *C, double *A, double *B, size_t dim);
void tensor_copy(double *C, double *A, size_t dim);
void tensor_decompose_deviatoric_and_spherical(double *A_0, double *c, double *A, size_t dim);

#endif // __TENSOR_H__
