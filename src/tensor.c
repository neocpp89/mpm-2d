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

#include "tensor.h"
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
void tensor_decompose(double *A_0, double *c, double *A, size_t dim)
{
    size_t i;
    double trA;
    assert(dim == 2 || dim == 3);
    assert(A != NULL);
    assert(A_0 != NULL);
    assert(c != NULL);
    assert(A != A_0);


    tensor_trace(&trA, A, dim);
    *c = trA / (double)dim;

    tensor_copy(A_0, A, dim);

    for (i = 0; i < dim; i++) {
        A_0[i*dim + i] = A_0[i*dim + i] - (*c);
    }

    return;
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \brief Multiples two tensors and stores the result in a third.
    \param C where to store result tensor.
    \param A left tensor.
    \param B right tensor.
    \param dim Is this a 2D or 3D tensor? Pass in either 2 or 3.
*/
void tensor_multiply(double *C, double *A, double *B, size_t dim)
{
    size_t i, j, k;
    assert(dim == 2 || dim == 3);
    assert(C != NULL);
    assert(A != NULL);
    assert(B != NULL);
    assert(A != C);

    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            C[i*dim + j] = 0;
            for (k = 0; k < dim; k++) {
                 C[i*dim + j] += A[i*dim + k]*B[k*dim + j];
            }
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \brief Scales tensor by a constant.
    \param A tensor to scale.
    \param c scaling constant.
    \param dim Is this a 2D or 3D tensor? Pass in either 2 or 3.
*/
void tensor_scale(double *A, double c, size_t dim)
{
    size_t i, j;
    assert(dim == 2 || dim == 3);
    assert(A != NULL);

    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            A[i*dim + j] *= c;
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \brief Gets symmetric part of tensor.
    \param C where to store sym(A).
    \param A tensor to symmetrize.
    \param dim Is this a 2D or 3D tensor? Pass in either 2 or 3.
*/
void tensor_sym(double *C, double *A, size_t dim)
{
    size_t i, j;
    assert(dim == 2 || dim == 3);
    assert(A != NULL);
    assert(C != NULL);
    assert(A != C);

    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            C[i*dim + j] = 0.5 *(A[i*dim + j] + A[j*dim + i]);
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \brief Gets skew part of tensor.
    \param C where to store skw(A).
    \param A tensor to skew.
    \param dim Is this a 2D or 3D tensor? Pass in either 2 or 3.
*/
void tensor_skw(double *C, double *A, size_t dim)
{
    size_t i, j;
    assert(dim == 2 || dim == 3);
    assert(A != NULL);
    assert(C != NULL);
    assert(A != C);

    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            C[i*dim + j] = 0.5 *(A[i*dim + j] - A[j*dim + i]);
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \brief Tensor contraction of A and B.
    \param C where to store skw(A).
    \param A tensor to skew.
    \param dim Is this a 2D or 3D tensor? Pass in either 2 or 3.
*/
void tensor_contraction(double *c, double *A, double *B, size_t dim)
{
    size_t i, j;
    assert(dim == 2 || dim == 3);
    assert(A != NULL);
    assert(B != NULL);
    assert(c != NULL);
    
    *c = 0;
    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            (*c) += A[i*dim + j] * B[i*dim + j];
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/
