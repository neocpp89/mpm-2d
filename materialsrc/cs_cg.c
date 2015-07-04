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
        s += (a[i] * b[i]);
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
    printf("n = %zu\n", n);

    double * restrict u = malloc(sizeof(double) * n);
    double * restrict r = malloc(sizeof(double) * n);
    double * restrict p = malloc(sizeof(double) * n);
    const double rsq_tol = tol*tol;

    if (u_0 != NULL) {
        memcpy(u, u_0, sizeof(double) * n);
        for (size_t j = 0; j < n; j++) {
            r[j] = 0;
        }

        cs_gaxpy(K, u, r); // r is now K * u_0
        for (size_t j = 0; j < n; j++) {
            r[j] = f[j] - r[j];
        }
    } else {
        for (size_t j = 0; j < n; j++) {
            u[j] = 0;
        }
        memcpy(r, f, sizeof(double) * n);
    }
    memcpy(p, r, sizeof(double) * n);

    /*
    for (size_t j = 0; j < n; j++) {
        printf("p[%zu] = %g\n", j, p[j]);
    }
    */

    double * restrict Kp = malloc(sizeof(double) * n);
    for (size_t i = 0; i < n*n*n; i++) {
        if (i % residual_recalculation_interval == 0) {
            for (size_t j = 0; j < n; j++) {
                r[j] = 0;
            }
            cs_gaxpy(K, u, r); // r is now K * u_i
            for (size_t j = 0; j < n; j++) {
                // printf("Ku[%zu] = %g\n", j, r[j]);
                r[j] = f[j] - r[j];
                // printf("f[%zu] = %g, ", j, f[j]);
                // printf("r[%zu] = %g\n", j, r[j]);
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
        /*
        if (pTKp == 0) {
            converged = (rTr == 0);
            break;
        }
        */
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

int cs_bicgstab(const cs *A, double *b, const double *x_0, double tol)
{
    const size_t residual_recalculation_interval = 10;
    const double tolsq = tol*tol;
    int converged = 0;
    if (A->m != A->n) {
        return converged; //nonsquare matrix
    }
    const size_t n = A->n;

    double *rhat_0 = calloc(n, sizeof(double));
    double *p_i = calloc(n, sizeof(double));
    double *v_i = calloc(n, sizeof(double));
    double *x_i = calloc(n, sizeof(double));
    double *r_i = calloc(n, sizeof(double));
    double *s = calloc(n, sizeof(double));
    double *t = calloc(n, sizeof(double));

    double bTb = dot(b, b, n);
    if (bTb == 0) {
        bTb = 1;
    } 

    if (x_0 != NULL) {
        cs_gaxpy(A, x_0, r_i); // r_i is Ax_0
        for (size_t j = 0; j < n; j++) {
            r_i[j] = b[j] - r_i[j]; // r_i is b - Ax_0
        }
    } else {
        for (size_t j = 0; j < n; j++) {
            r_i[j] = b[j];
        }
    }

    for (size_t j = 0; j < n; j++) {
        rhat_0[j] = r_i[j];
    }

    double error = dot(r_i, r_i, n) / bTb;
    if (error < tolsq) {
        converged = 1;
        goto freetemp;
    }

    double rho_i = 1;
    double omega_i = 1;
    double alpha = 1; 

    for (size_t j = 0; j < n; j++) {
        v_i[j] = 0;
        p_i[j] = 0;
    }

    for (size_t i = 0; i < n; i++) {
        const double rho_prev = rho_i;
        rho_i = dot(rhat_0, r_i, n);
        const double omega_prev = omega_i;
        const double beta = (rho_i / rho_prev) * (alpha / omega_prev);
        for (size_t j = 0; j < n; j++) {
            p_i[j] = r_i[j] + beta * (p_i[j] - omega_prev * v_i[j]);
        }
        for (size_t j = 0; j < n; j++) {
            v_i[j] = 0;
        }
        cs_gaxpy(A, p_i, v_i);
        alpha = rho_i / dot(rhat_0, v_i, n);
        for (size_t j = 0; j < n; j++) {
            s[j] = r_i[j] - alpha * v_i[j];
        }
        for (size_t j = 0; j < n; j++) {
            t[j] = 0;
        }
        cs_gaxpy(A, s, t);
        omega_i = dot(s, t, n) / dot(t, t, n);
        for (size_t j = 0; j < n; j++) {
            x_i[j] = x_i[j] + alpha * p_i[j] + omega_i * s[j];
        }
        for (size_t j = 0; j < n; j++) {
            r_i[j] = s[j] - omega_i * t[j];
        }
        const double rtr = dot(r_i, r_i, n);
        if (rtr < tolsq) {
            converged = 1;
            break;
        }
    }

freetemp:

    free(rhat_0);
    free(p_i);
    free(v_i);
    free(x_i);
    free(r_i);
    free(s);
    free(t);
    return converged;
}
