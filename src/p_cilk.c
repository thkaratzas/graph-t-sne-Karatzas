#include <stdlib.h>
#include <cilk/cilk.h>
#include "tsne.h"

#define GRAIN 256

void compute_P_cilk(const graph_t *W, double u, double *beta, double *P)
{
    int64_t n = W->n;
    double *pc = malloc(W->nnz * sizeof(double));   /* p_{j|i} */

    /* conditional probabilities, row by row.
       The row's own part of pc is also its work buffer (see row_beta). */
    #pragma cilk grainsize GRAIN
    cilk_for (int64_t i = 0; i < n; i++) {
        int64_t k0 = W->colptr[i];
        int64_t d  = W->colptr[i + 1] - k0;
        beta[i] = row_beta(&W->val[k0], d, u, &pc[k0], &pc[k0]);
    }

    /* symmetrize: p_ij = (p_{j|i} + p_{i|j}) / 2n */
    #pragma cilk grainsize GRAIN
    cilk_for (int64_t i = 0; i < n; i++) {
        for (int64_t k = W->colptr[i]; k < W->colptr[i + 1]; k++) {
            int64_t j  = W->row[k];
            int64_t kk = find_entry(W, j, i);
            P[k] = (pc[k] + pc[kk]) / (2.0 * n);
        }
    }

    free(pc);
}