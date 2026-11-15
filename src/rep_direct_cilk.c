#include <cilk/cilk.h>
#include "tsne.h"

static void zero_d(void *v)          { *(double *)v = 0.0; }
static void plus_d(void *l, void *r) { *(double *)l += *(double *)r; }

/* Direct O(n^2) repulsion, parallel over i. Returns Z. */
double compute_rep_direct_cilk(int64_t n, const double *Y, double *F)
{
    double cilk_reducer(zero_d, plus_d) Z = 0;

    cilk_for (int64_t i = 0; i < n; i++) {
        double yi0 = Y[2 * i];
        double yi1 = Y[2 * i + 1];
        double f0 = 0, f1 = 0, zi = 0;

        for (int64_t j = 0; j < n; j++) {
            if (j == i)
                continue;
            double dx = yi0 - Y[2 * j];
            double dy = yi1 - Y[2 * j + 1];
            double q  = 1.0 / (1.0 + dx * dx + dy * dy);
            zi += q;
            f0 += q * q * dx;
            f1 += q * q * dy;
        }

        F[2 * i]     = f0;
        F[2 * i + 1] = f1;
        Z += zi;
    }

    double Ztot = Z;

    cilk_for (int64_t k = 0; k < 2 * n; k++)
        F[k] /= Ztot;

    return Ztot;
}