#include "tsne.h"

/* F_attr_i = sum_j p_ij (1 + |y_i - y_j|^2)^-1 (y_i - y_j) */
void compute_attr_seq(const graph_t *W, const double *P, const double *Y, double *F)
{
    for (int64_t i = 0; i < W->n; i++) {
        double yi0 = Y[2 * i];
        double yi1 = Y[2 * i + 1];
        double f0 = 0, f1 = 0;

        for (int64_t k = W->colptr[i]; k < W->colptr[i + 1]; k++) {
            int64_t j = W->row[k];
            double dx = yi0 - Y[2 * j];
            double dy = yi1 - Y[2 * j + 1];
            double a  = P[k] / (1.0 + dx * dx + dy * dy);
            f0 += a * dx;
            f1 += a * dy;
        }

        F[2 * i]     = f0;
        F[2 * i + 1] = f1;
    }
}