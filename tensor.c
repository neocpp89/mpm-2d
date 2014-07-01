/**
    \file tensor.c
    \author Sachith Dunatunga
    \date 04.06.12

    mpm_2d -- An implementation of the Material Point Method in 2D.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "particle.h"

#include <assert.h>

/*----------------------------------------------------------------------------*/
void zero_array(double *array, size_t numel)
{
#ifdef __STDC_IEC_559__
    /*
        We conform to IEEE754-1985, so we can safely use memset to set all
        zero bytes (0b -> 0.0).
    */
    memset(array, 0, numel * sizeof(double));
#else
    size_t i;
    for (i = 0; i < numel; i++) {
        array[i] = 0;
    }
#endif
    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/** \brief Calculates the deviatoric and spherical part of a tensor.
    \param A_0 Where to store the deviator.
    \param A Tensor to take the deviator of.
    \param dim Is this a 2D or 3D tensor? Pass in either 2 or 3.
*/
void tensor_trace(double *trA, double *A, size_t dim)
{
    size_t i;
    assert(dim == 2 || dim == 3);

/*    if (dim == 2) {*/
/*        *trA = A[XX] + A[YY];*/
/*    } else if (dim == 3) {*/
/*        *trA = A[XX] + A[YY] + A[ZZ];*/
/*    } else {*/
/*        *trA = 0;*/
/*    }*/
    *trA = 0;
    for (i = 0; i < dim; i++) {
        *trA += A[i * dim + i];
    }
    
    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/** \brief Calculates the deviatoric and spherical part of a tensor.
    \param C Where to store the result (tensor).
    \param A First tensor (tensor).
    \param B Second tensor (tensor).
    \param dim Is this a 2D or 3D tensor? Pass in either 2 or 3.

    This function adds the B tensor to the A tensor and stores the result in
    a tensor C. You may use the same tensor for both A and C to overwrite the
    value of A.
*/
void tensor_add(double *C, double *A, double *B, size_t dim)
{
    size_t i, j;
    assert(dim == 2 || dim == 3);
    assert(C != NULL);
    assert(A != NULL);
    assert(B != NULL);

    if (C != A) {
        tensor_copy(C, A, dim*dim);
    }

    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            C[i*dim + j] += B[i*dim + j];
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/** \brief Calculates the deviatoric and spherical part of a tensor.
    \param C Where to store the result (tensor).
    \param A Tensor to copy (tensor).
    \param dim Is this a 2D or 3D tensor? Pass in either 2 or 3.
*/
void tensor_copy(double *C, double *A, size_t dim)
{
    size_t i, j;
    assert(dim == 2 || dim == 3);
    assert(C != NULL);
    assert(A != NULL);

    if (C == A) {
        return;
    }

    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            C[i*dim + j] = A[i*dim + j];
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/** \brief Calculates the deviatoric and spherical part of a tensor.
    \param A_0 Where to store the deviator (tensor).
    \param cI Where to store the spherical part of A (scalar).
    \param A Tensor to take the deviator of.
    \param dim Is this a 2D or 3D tensor? Pass in either 2 or 3.
*/
void tensor_decompose_deviatoric_and_spherical(double *A_0, double *c, double *A, size_t dim)
{
    size_t i;
    double trA;
    assert(dim == 2 || dim == 3);
    assert(A != NULL);
    assert(A_0 != NULL);
    assert(c != NULL);
    assert(A != A_0);


    tensor_trace(&trA, A, dim);
    *c = trA / (float)dim;

    zero_array(A_0, dim*dim);
    tensor_copy(A_0, A, dim);

    for (i = 0; i < dim; i++) {
        A_0[i*dim + i] = A_0[i*dim + i] - (*c);
    }

    return;
}
/*----------------------------------------------------------------------------*/

