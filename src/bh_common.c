#include "tsne.h"

/* Barnes-Hut walk for one point i.
   acc[0], acc[1]: unnormalized repulsion of point i,  acc[2]: its part of Z */
void bh_walk(const cell_t *c, const int64_t *perm, const double *Y,
             int64_t i, double theta, double *acc)
{
    double yi0 = Y[2 * i];
    double yi1 = Y[2 * i + 1];

    /* far enough: the whole cell as one point at its center of mass */
    double dx = yi0 - c->mx;
    double dy = yi1 - c->my;
    double d2 = dx * dx + dy * dy;
    double r  = 2.0 * c->h * 1.4142135623730951;   /* diagonal of the cell */
    if (r * r < theta * theta * d2) {
        double q = 1.0 / (1.0 + d2);
        acc[0] += c->count * q * q * dx;
        acc[1] += c->count * q * q * dy;
        acc[2] += c->count * q;
        return;
    }

    /* leaf: every point, except i itself */
    if (c->leaf) {
        for (int64_t k = c->start; k < c->start + c->count; k++) {
            int64_t j = perm[k];
            if (j == i)
                continue;
            double ex = yi0 - Y[2 * j];
            double ey = yi1 - Y[2 * j + 1];
            double q  = 1.0 / (1.0 + ex * ex + ey * ey);
            acc[0] += q * q * ex;
            acc[1] += q * q * ey;
            acc[2] += q;
        }
        return;
    }

    /* open the cell */
    for (int q = 0; q < 4; q++)
        if (c->child[q])
            bh_walk(c->child[q], perm, Y, i, theta, acc);
}