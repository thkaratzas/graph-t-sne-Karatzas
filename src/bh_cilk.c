#include <cilk/cilk.h>
#include "tsne.h"

#define GRAIN 256

/* reducer for a sum of doubles:
   identity sets a new view to 0, reduce adds the right view into the left */
static void zero_d(void *v)          { *(double *)v = 0.0; }
static void plus_d(void *l, void *r) { *(double *)l += *(double *)r; }

double compute_rep_bh_cilk(const tree_t *T, int64_t n, const double *Y, double theta, double *F)
{
    double cilk_reducer(zero_d, plus_d) Z = 0;

    #pragma cilk grainsize GRAIN
    cilk_for (int64_t i = 0; i < n; i++) {
        double acc[3] = {0, 0, 0};
        bh_walk(T->root, T->perm, Y, i, theta, acc);
        F[2 * i]     = acc[0];
        F[2 * i + 1] = acc[1];
        Z += acc[2];
    }

    double Ztot = Z;

    cilk_for (int64_t k = 0; k < 2 * n; k++)
        F[k] /= Ztot;

    return Ztot;
}