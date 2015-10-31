#ifndef SPMD_H
#define SPMD_H

struct sparsematrix_double {
    double *vals;
    size_t nnz;

    size_t *column_index;
    size_t *row_pointer;
    size_t rows;
    size_t columns;
};

struct sparsematrix_double *spmd_create(size_t nnz, size_t rows, size_t columns);
void spmd_delete(struct sparsematrix_double *sp);
void spmd_print(const struct sparsematrix_double *sp, int full);

void spmd_gaxpy(const struct sparsematrix_double *A, const double *x, double *y);
void spmdv(double *result, const struct sparsematrix_double *sp, const double *v);

#endif // SPMD_H
