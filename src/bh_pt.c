#include <stdlib.h>
#include "tsne.h"

/* one partial sum per thread, padded to a full cache line (64 bytes)
   so that two threads never write to the same line */
typedef struct {
    double z;
    char pad[64 - sizeof(double)];
} partial_t;

typedef struct {
    const tree_t *T;
    const double *Y;
    double theta;
    double *F;
    partial_t *partial;   /* [number of threads] */
    double Z;             /* set before the normalization loop */
} bh_args_t;

/* walk for points lo..hi-1 */
static void walk_body(int64_t lo, int64_t hi, int tid, void *v)
{
    bh_args_t *a = v;
    double z = 0;
    for (int64_t i = lo; i < hi; i++) {
        double acc[3] = {0, 0, 0};
        bh_walk(a->T->root, a->T->perm, a->Y, i, a->theta, acc);
        a->F[2 * i]     = acc[0];
        a->F[2 * i + 1] = acc[1];
        z += acc[2];
    }
    a->partial[tid].z += z;   /* only this thread writes here */
}

/* F_k /= Z for k in lo..hi-1 */
static void norm_body(int64_t lo, int64_t hi, int tid, void *v)
{
    bh_args_t *a = v;
    for (int64_t k = lo; k < hi; k++)
        a->F[k] /= a->Z;
}

double compute_rep_bh_pt(const tree_t *T, int64_t n, const double *Y, double theta, double *F)
{
    int p = pt_nthreads();
    bh_args_t a;
    a.T = T;
    a.Y = Y;
    a.theta = theta;
    a.F = F;
    a.partial = calloc(p, sizeof(partial_t));   /* all zero */

    pt_for(n, walk_body, &a);

    double Z = 0;
    for (int t = 0; t < p; t++)
        Z += a.partial[t].z;
    a.Z = Z;

    pt_for(2 * n, norm_body, &a);

    free(a.partial);
    return Z;
}