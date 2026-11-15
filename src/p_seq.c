#include <stdlib.h>
#include "tsne.h"

void compute_P_seq(const graph_t *W, double u, double *beta, double *P)
{
    int64_t n = W->n;
    double *pc = malloc(W->nnz * sizeof(double));          /* p_{j|i} */
    double *lr = malloc(max_degree(W) * sizeof(double));   /* work buffer */

    /* conditional probabilities, row by row */
    for (int64_t i = 0; i < n; i++) {
        int64_t k0 = W->colptr[i];
        int64_t d  = W->colptr[i + 1] - k0;
        beta[i] = row_beta(&W->val[k0], d, u, &pc[k0], lr);
    }

    /* symmetrize: p_ij = (p_{j|i} + p_{i|j}) / 2n */
    for (int64_t i = 0; i < n; i++) {
        for (int64_t k = W->colptr[i]; k < W->colptr[i + 1]; k++) {
            int64_t j  = W->row[k];
            int64_t kk = find_entry(W, j, i);
            P[k] = (pc[k] + pc[kk]) / (2.0 * n);
        }
    }

    free(pc);
    free(lr);
}