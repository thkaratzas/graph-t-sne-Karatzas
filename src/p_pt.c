#include <stdlib.h>
#include "tsne.h"

/* data that both loops need */
typedef struct {
    const graph_t *W;
    double u;
    double *beta;
    double *pc;    /* p_{j|i} */
    double *P;
} p_args_t;

/* rows lo..hi-1: conditional probabilities.
   The row's own part of pc is also its work buffer (as in p_cilk.c). */
static void rows_body(int64_t lo, int64_t hi, int tid, void *v)
{
    p_args_t *a = v;
    for (int64_t i = lo; i < hi; i++) {
        int64_t k0 = a->W->colptr[i];
        int64_t d  = a->W->colptr[i + 1] - k0;
        a->beta[i] = row_beta(&a->W->val[k0], d, a->u, &a->pc[k0], &a->pc[k0]);
    }
}

/* rows lo..hi-1: p_ij = (p_{j|i} + p_{i|j}) / 2n */
static void sym_body(int64_t lo, int64_t hi, int tid, void *v)
{
    p_args_t *a = v;
    const graph_t *W = a->W;
    for (int64_t i = lo; i < hi; i++) {
        for (int64_t k = W->colptr[i]; k < W->colptr[i + 1]; k++) {
            int64_t j  = W->row[k];
            int64_t kk = find_entry(W, j, i);
            a->P[k] = (a->pc[k] + a->pc[kk]) / (2.0 * W->n);
        }
    }
}

void compute_P_pt(const graph_t *W, double u, double *beta, double *P)
{
    p_args_t a;
    a.W = W;
    a.u = u;
    a.beta = beta;
    a.pc = malloc(W->nnz * sizeof(double));
    a.P = P;

    pt_for(W->n, rows_body, &a);   /* all rows done when it returns */
    pt_for(W->n, sym_body, &a);

    free(a.pc);
}