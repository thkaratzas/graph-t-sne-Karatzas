#include <stdlib.h>
#include "tsne.h"

typedef struct {
    double z;
    char pad[64 - sizeof(double)];
} partial_t;

typedef struct {
    int64_t n;
    const double *Y;
    double *F;
    partial_t *partial;
    double Z;
} dir_args_t;

/* direct repulsion for points lo..hi-1 */
static void dir_body(int64_t lo, int64_t hi, int tid, void *v)
{
    dir_args_t *a = v;
    const double *Y = a->Y;
    double z = 0;

    for (int64_t i = lo; i < hi; i++) {
        double yi0 = Y[2 * i];
        double yi1 = Y[2 * i + 1];
        double f0 = 0, f1 = 0, zi = 0;

        for (int64_t j = 0; j < a->n; j++) {
            if (j == i)
                continue;
            double dx = yi0 - Y[2 * j];
            double dy = yi1 - Y[2 * j + 1];
            double q  = 1.0 / (1.0 + dx * dx + dy * dy);
            zi += q;
            f0 += q * q * dx;
            f1 += q * q * dy;
        }

        a->F[2 * i]     = f0;
        a->F[2 * i + 1] = f1;
        z += zi;
    }
    a->partial[tid].z += z;
}

static void norm_body(int64_t lo, int64_t hi, int tid, void *v)
{
    dir_args_t *a = v;
    for (int64_t k = lo; k < hi; k++)
        a->F[k] /= a->Z;
}

double compute_rep_direct_pt(int64_t n, const double *Y, double *F)
{
    int p = pt_nthreads();
    dir_args_t a;
    a.n = n;
    a.Y = Y;
    a.F = F;
    a.partial = calloc(p, sizeof(partial_t));

    pt_for(n, dir_body, &a);

    double Z = 0;
    for (int t = 0; t < p; t++)
        Z += a.partial[t].z;
    a.Z = Z;

    pt_for(2 * n, norm_body, &a);

    free(a.partial);
    return Z;
}