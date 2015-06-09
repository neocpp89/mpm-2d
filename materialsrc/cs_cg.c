/**
    \file cs_cg.c
    \author Sachith Dunatunga
    \date 06.05.2015

    mpm_2d -- An implementation of the Material Point Method in 2D.
*/
#include <suitesparse/cs.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "cs_cg.h"

double dot(const double * a, const double * b, size_t n)
{
    double s = 0;
    for (size_t i = 0; i < n; i++) {
        s += a[i] * b[i];
    }
    return s;
}

int cs_cg(const cs *K, double *f, const double *u_0, double tol)
{
    const size_t residual_recalculation_interval = 10;
    int converged = 0;
    if (K->m != K->n) {
        return converged; //nonsquare matrix
    }
    const size_t n = K->n;

    double * restrict u = malloc(sizeof(double) * n);
    double * restrict r = malloc(sizeof(double) * n);
    double * restrict p = malloc(sizeof(double) * n);
    const double rsq_tol = tol*tol;

    if (u_0 != NULL) {
        memcpy(u, u_0, sizeof(double) * n);
        for (size_t j = 0; j < n; j++) {
            r[j] = 0;
        }
        // Taken care of once we start the iteration.
        /* 
        cs_gaxpy(K, u, r); // r is now K * u_0
        for (size_t j = 0; j < n; j++) {
            r[j] = f[j] - r[j];
        }
        */
    } else {
        for (size_t j = 0; j < n; j++) {
            u[j] = 0;
        }
        memcpy(r, f, sizeof(double) * n);
    }
    memcpy(p, r, sizeof(double) * n);

    double * restrict Kp = malloc(sizeof(double) * n);
    for (size_t i = 0; i < n; i++) {
        if (i % residual_recalculation_interval == 0) {
            cs_gaxpy(K, u, r); // r is now K * u_i
            for (size_t j = 0; j < n; j++) {
                r[j] = f[j] - r[j];
            }
        }

        const double rTr = dot(r, r, n);
        // printf("r.r = %lg\n", rTr);
        if (rTr <= rsq_tol) {
            converged = 1;
            // printf("\n");
            break;
        }
        for (size_t j = 0; j < n; j++) {
            Kp[j] = 0;
        }
        cs_gaxpy(K, p, Kp); // Kp is now K * p_{i-1}
        const double pTKp = dot(p, Kp, n);
        assert(pTKp != 0);
        const double alpha = rTr / pTKp;
        for (size_t j = 0; j < n; j++) {
            u[j] = u[j] + alpha * p[j];
            r[j] = r[j] - alpha * Kp[j];
        }
        const double rTr_new = dot(r, r, n);
        const double beta = rTr_new / rTr;
        for (size_t j = 0; j < n; j++) {
            p[j] = r[j] + beta * p[j];
        }
    }

    memcpy(f, u, sizeof(double) * n);
    free(u);
    free(r);
    free(p);
    free(Kp);
    return converged;
}

