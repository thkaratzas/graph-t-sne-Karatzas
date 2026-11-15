#include <omp.h>
#include "tsne.h"

/* Step 4: grad_i = 4 (Fattr_i - Frep_i) */
void compute_grad_omp(int64_t n, const double *Fattr, const double *Frep, double *G)
{
    #pragma omp parallel for schedule(static)
    for (int64_t k = 0; k < 2 * n; k++)
        G[k] = 4.0 * (Fattr[k] - Frep[k]);
}