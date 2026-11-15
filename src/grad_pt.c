#include "tsne.h"

typedef struct {
    const double *Fattr;
    const double *Frep;
    double *G;
} grad_args_t;

/* G_k = 4 (Fattr_k - Frep_k), for k in lo..hi-1 (k runs over 2n coordinates) */
static void grad_body(int64_t lo, int64_t hi, int tid, void *v)
{
    grad_args_t *a = v;
    for (int64_t k = lo; k < hi; k++)
        a->G[k] = 4.0 * (a->Fattr[k] - a->Frep[k]);
}

void compute_grad_pt(int64_t n, const double *Fattr, const double *Frep, double *G)
{
    grad_args_t a = {Fattr, Frep, G};
    pt_for(2 * n, grad_body, &a);
}