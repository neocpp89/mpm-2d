#include <stdio.h>
#include <stdlib.h>

#include "spmd.h"

void do_netlib_example();
void do_rectangular_example();
void do_repeated_index();

void do_netlib_example()
{
    struct sparsematrix_double *sp = spmd_create(19, 6, 6);
    const double vals[19] = {10, -2, 3, 9, 3, 7, 8, 7, 3, 8, 7, 5, 8, 9, 9, 13, 4, 2, -1};
    const size_t column_index[19] = {0, 4, 0, 1, 5, 1, 2, 3, 0, 2, 3, 4, 1, 3, 4, 5, 1, 4, 5};
    const size_t row_pointer[6] = {0, 2, 5, 8, 12, 16};
    for (size_t i = 0; i < sp->nnz; i++) {
        sp->vals[i] = vals[i];
        sp->column_index[i] = column_index[i];
    }
    for (size_t i = 0; i < sp->rows; i++) {
        sp->row_pointer[i] = row_pointer[i];
    }
    spmd_print(sp, 1);

    double result[6] = {0, 0, 0, 0, 0, 0};
    double x[6] = {1, 1, 1, 1, 1, 1};
    spmd_gaxpy(sp, x, result);
    printf("x: ");
    for (size_t i = 0; i < sp->rows; i++) {
        printf("%4.4g ", x[i]);
    }
    printf("\n");
    printf("result: ");
    for (size_t i = 0; i < sp->rows; i++) {
        printf("%4.4g ", result[i]);
    }
    printf("\n");

    double tresult[6] = {0, 0, 0, 0, 0, 0};
    spmd_gatxpy(sp, x, tresult);
    printf("tresult: ");
    for (size_t i = 0; i < sp->columns; i++) {
        printf("%4.4g ", tresult[i]);
    }
    printf("\n");

    spmd_delete(sp);
}

void do_rectangular_example()
{
    struct sparsematrix_double *sp = spmd_create(5, 2, 4);
    const double vals[] = {4, -7, 1, 1, 2};
    const double column_index[] = {1, 2, 3, 0, 3};
    const double row_pointer[] = {0, 3};
    for (size_t i = 0; i < sp->nnz; i++) {
        sp->vals[i] = vals[i];
        sp->column_index[i] = column_index[i];
    }
    for (size_t i = 0; i < sp->rows; i++) {
        sp->row_pointer[i] = row_pointer[i];
    }
    spmd_print(sp, 1);

    double result[] = {0, 0};
    double x[] = {1, 1, 1, 1};
    spmd_gaxpy(sp, x, result);
    printf("x: ");
    for (size_t i = 0; i < sp->rows; i++) {
        printf("%4.4g ", x[i]);
    }
    printf("\n");
    printf("result: ");
    for (size_t i = 0; i < sp->rows; i++) {
        printf("%4.4g ", result[i]);
    }
    printf("\n");

    double tresult[] = {0, 0, 0, 0};
    double tx[] = {1, 1};
    spmd_gatxpy(sp, tx, tresult);
    printf("tresult: ");
    for (size_t i = 0; i < sp->columns; i++) {
        printf("%4.4g ", tresult[i]);
    }
    printf("\n");

    spmd_delete(sp);
}

void do_repeated_index()
{
    struct sparsematrix_double *sp = spmd_create(6, 2, 4);
    const double vals[] = {4, -7, 1, 1, 2, 10};
    const double column_index[] = {1, 2, 3, 0, 3, 3};
    const double row_pointer[] = {0, 3};
    for (size_t i = 0; i < sp->nnz; i++) {
        sp->vals[i] = vals[i];
        sp->column_index[i] = column_index[i];
    }
    for (size_t i = 0; i < sp->rows; i++) {
        sp->row_pointer[i] = row_pointer[i];
    }
    spmd_print(sp, 1);

    double result[] = {0, 0};
    double x[] = {1, 1, 1, 1};
    spmd_gaxpy(sp, x, result);
    printf("x: ");
    for (size_t i = 0; i < sp->rows; i++) {
        printf("%4.4g ", x[i]);
    }
    printf("\n");
    printf("result: ");
    for (size_t i = 0; i < sp->rows; i++) {
        printf("%4.4g ", result[i]);
    }
    printf("\n");

    double tresult[] = {0, 0, 0, 0};
    double tx[] = {1, 1};
    spmd_gatxpy(sp, tx, tresult);
    printf("tresult: ");
    for (size_t i = 0; i < sp->columns; i++) {
        printf("%4.4g ", tresult[i]);
    }
    printf("\n");

    spmd_delete(sp);
}
int main(int argc, char **argv)
{
    do_netlib_example();
    do_rectangular_example();
    do_repeated_index();

    return 0;
}
