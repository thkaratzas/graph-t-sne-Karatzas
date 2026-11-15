#include <math.h>
#include "tsne.h"

#define BETA_TOL   1e-5
#define BETA_MAXIT 200

/* Entropy of p_k = exp(beta*lr_k) / S */
static double entropy(const double *lr, int64_t d, double beta)
{
    double S = 0, T = 0;
    for (int64_t k = 0; k < d; k++) {
        double e = exp(beta * lr[k]);
        S += e;
        T += e * lr[k];
    }
    return log(S) - beta * T / S;
}

/* One row: weights w[0..d-1] > 0.
   Writes p[k] = p_{j|i} and returns beta. lr is a work buffer of size >= d. */
double row_beta(const double *w, int64_t d, double u, double *p, double *lr)
{
    double target = log(d < u ? (double)d : u);

    double wmax = w[0];
    for (int64_t k = 1; k < d; k++)
        if (w[k] > wmax)
            wmax = w[k];

    int all_equal = 1;
    for (int64_t k = 0; k < d; k++) {
        lr[k] = log(w[k] / wmax);
        if (lr[k] != 0.0)
            all_equal = 0;
    }

    double beta = 0;
    if (!all_equal && u < d) {
        double lo = 0, hi = -1;   /* hi < 0 means: no upper bound yet */
        beta = 1;
        for (int it = 0; it < BETA_MAXIT; it++) {
            double H = entropy(lr, d, beta);
            if (fabs(H - target) < BETA_TOL)
                break;
            if (H > target) {            /* too flat: increase beta */
                lo = beta;
                beta = (hi < 0) ? 2 * beta : (lo + hi) / 2;
            } else {                     /* too sharp: decrease beta */
                hi = beta;
                beta = (lo + hi) / 2;
            }
        }
    }

    double S = 0;
    for (int64_t k = 0; k < d; k++) {
        p[k] = exp(beta * lr[k]);
        S += p[k];
    }
    for (int64_t k = 0; k < d; k++)
        p[k] /= S;

    return beta;
}

/* Position of entry (j, i) in row j, or -1. Rows are sorted. */
int64_t find_entry(const graph_t *W, int64_t j, int64_t i)
{
    int64_t lo = W->colptr[j];
    int64_t hi = W->colptr[j + 1] - 1;
    while (lo <= hi) {
        int64_t mid = (lo + hi) / 2;
        if (W->row[mid] == i)
            return mid;
        if (W->row[mid] < i)
            lo = mid + 1;
        else
            hi = mid - 1;
    }
    return -1;
}

int64_t max_degree(const graph_t *W)
{
    int64_t dmax = 0;
    for (int64_t i = 0; i < W->n; i++) {
        int64_t d = W->colptr[i + 1] - W->colptr[i];
        if (d > dmax)
            dmax = d;
    }
    return dmax;
}