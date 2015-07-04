#include <stdio.h>
#include <stdlib.h>

#include "cs_cg.h"

void poisson_matrix(cs *triplets, size_t lda)
{
    const double h = 1.0 / (lda-1);
    const double hinv_sq = 1.0 / (h * h);
    cs_entry(triplets, 0, 0, 1);
    for (size_t i = 1; i < (lda-1); i++) {
        cs_entry(triplets, i, i-1, 1 * hinv_sq);
        cs_entry(triplets, i, i, -2* hinv_sq);
        cs_entry(triplets, i, i+1, 1 * hinv_sq);
    }
    cs_entry(triplets, lda-1, lda-1, 1);
    return;
}

int main(int argc, char **argv)
{
    const size_t lda = 101;
    const size_t nnz = 3 * (lda - 2) + 2;
    cs *triplets = cs_spalloc(lda, lda, nnz, 1, 1);
    printf("n=%d, m=%d\n", triplets->n, triplets->m);
    poisson_matrix(triplets, lda);
    cs *A = cs_compress(triplets);
    cs_dupl(A);
    printf("n=%d, m=%d\n", A->n, A->m);
    double *b = calloc(lda, sizeof(double));
    double *z = calloc(lda, sizeof(double));
    b[lda-1] = 1;

    for (size_t i = 0; i < lda; i++) {
        // z[i] = (double)i / (double)(lda-1);
    }

    int r = cs_cg(A, b, b, 1e-10);
    // int r = cs_bicgstab(A, b, b, 1e-10);
    // int r = cs_lusol(1, A, b, 1e-10);

    for (size_t i = 0; i < lda; i++) {
        printf("x[%zu] = %g\n", i, b[i]);
    }

    cs_spfree(triplets);
    cs_spfree(A);
    free(b);
    free(z);
    if (r == 0) {
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
