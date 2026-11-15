#include "tsne.h"

double compute_rep_bh_seq(const tree_t *T, int64_t n, const double *Y, double theta, double *F)
{
    double Z = 0;

    for (int64_t i = 0; i < n; i++) {
        double acc[3] = {0, 0, 0};
        bh_walk(T->root, T->perm, Y, i, theta, acc);
        F[2 * i]     = acc[0];
        F[2 * i + 1] = acc[1];
        Z += acc[2];
    }

    for (int64_t k = 0; k < 2 * n; k++)
        F[k] /= Z;

    return Z;
}