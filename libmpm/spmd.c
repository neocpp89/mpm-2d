#include <stdlib.h>
#include <stdio.h>

#include "spmd.h"

double spmd_slow_get(const struct sparsematrix_double *sp, size_t i, size_t j);

struct sparsematrix_double *spmd_create(size_t nnz, size_t rows, size_t columns)
{
    struct sparsematrix_double *sp = malloc(sizeof(struct sparsematrix_double));
    if (sp == NULL) {
        return NULL;
    }
    sp->nnz = nnz;
    sp->capacity = nnz;
    sp->rows = rows;
    sp->columns = columns;
    sp->vals = calloc(nnz, sizeof(double));
    if (sp->vals == NULL) {
        free(sp);
        return NULL;
    }
    sp->column_index = calloc(nnz, sizeof(size_t));
    if (sp->column_index == NULL) {
        free(sp->vals);
        free(sp);
        return NULL;
    }

    sp->row_pointer = calloc(rows + 1, sizeof(size_t));
    if (sp->row_pointer == NULL) {
        free(sp->column_index);
        free(sp->vals);
        free(sp);
        return NULL;
    }
    sp->row_pointer[rows] = nnz;

    return sp;
}

void spmd_delete(struct sparsematrix_double *sp)
{
    if (sp == NULL) {
        return;
    }

    free(sp->row_pointer);
    free(sp->column_index);
    free(sp->vals);
    free(sp);
    return;
}

void spmd_print(const struct sparsematrix_double *sp, int full)
{
    if (sp == NULL) {
        return;
    }

    printf("nnz: %zu, rows: %zu, cols: %zu\n", sp->nnz, sp->rows, sp->columns);
    if (full == 0) {
        printf("vals: ");
        for (size_t i = 0; i < sp->nnz; i++) {
            printf("%g ", sp->vals[i]);
        }
        printf("\ncolumn_index: ");
        for (size_t i = 0; i < sp->nnz; i++) {
            printf("%zu ", sp->column_index[i]);
        }
        printf("\nrow_pointer: ");
        for (size_t i = 0; i < sp->rows; i++) {
            printf("%zu ", sp->row_pointer[i]);
        }
        printf("\n");
    } else {
        for (size_t i = 0; i < sp->rows; i++) {
            for (size_t j = 0; j < sp->columns; j++) {
                printf("%4.4g ", spmd_slow_get(sp, i, j));
            }
            printf("\n");
        }
    }

    return;
}

double spmd_slow_get(const struct sparsematrix_double *sp, size_t i, size_t j)
{
    if (sp == NULL) {
        return 0xdeadbeef;
    }

    for (size_t p = sp->row_pointer[i]; p < sp->row_pointer[i+1]; p++) {
        if (sp->column_index[p] == j) {
            return sp->vals[p];
        }
    }

    return 0;
}

void spmd_gaxpy(const struct sparsematrix_double *A, const double *x, double *y)
{
    // does y <- A * x + y
    if (A == NULL || x == NULL || y == NULL) {
        return;
    }

    for (size_t i = 0; i < A->rows; i++) {
        for (size_t p = A->row_pointer[i]; p < A->row_pointer[i+1]; p++) {
            size_t j = A->column_index[p];
            y[i] += A->vals[p] * x[j];
        }
    }

    return;
}

void spmd_gatxpy(const struct sparsematrix_double *A, const double *x, double *y)
{
    // does y <- A^T * x + y
    if (A == NULL || x == NULL || y == NULL) {
        return;
    }

    for (size_t i = 0; i < A->rows; i++) {
        for (size_t p = A->row_pointer[i]; p < A->row_pointer[i+1]; p++) {
            size_t j = A->column_index[p];
            y[j] += A->vals[p] * x[i];
        }
    }

    return;
}

void spmdv(double *result, const struct sparsematrix_double *sp, const double *v)
{
    if (result == NULL || sp == NULL || v == NULL) {
        return;
    }

    for (size_t i = 0; i < sp->rows; i++) {
        result[i] = 0;
    }
    spmd_gaxpy(sp, v, result);

    return;
}
