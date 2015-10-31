#include <stdio.h>
#include <stdlib.h>

#include "spmd.h"

int main(int argc, char **argv)
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

    return 0;
}
