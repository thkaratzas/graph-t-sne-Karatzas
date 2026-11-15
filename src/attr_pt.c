#include "tsne.h"

typedef struct {
    const graph_t *W;
    const double *P;
    const double *Y;
    double *F;
} attr_args_t;

/* F_attr_i = sum_j p_ij (1 + |y_i - y_j|^2)^-1 (y_i - y_j), for i in lo..hi-1 */
static void attr_body(int64_t lo, int64_t hi, int tid, void *v)
{
    attr_args_t *a = v;
    const graph_t *W = a->W;
    const double *Y = a->Y;

    for (int64_t i = lo; i < hi; i++) {
        double yi0 = Y[2 * i];
        double yi1 = Y[2 * i + 1];
        double f0 = 0, f1 = 0;

        for (int64_t k = W->colptr[i]; k < W->colptr[i + 1]; k++) {
            int64_t j = W->row[k];
            double dx = yi0 - Y[2 * j];
            double dy = yi1 - Y[2 * j + 1];
            double w  = a->P[k] / (1.0 + dx * dx + dy * dy);
            f0 += w * dx;
            f1 += w * dy;
        }

        a->F[2 * i]     = f0;
        a->F[2 * i + 1] = f1;
    }
}

void compute_attr_pt(const graph_t *W, const double *P, const double *Y, double *F)
{
    attr_args_t a = {W, P, Y, F};
    pt_for(W->n, attr_body, &a);
}