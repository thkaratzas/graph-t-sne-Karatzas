#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "tsne.h"

static double now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}

/* sorted rows, symmetric pattern, positive weights */
static int check_graph(const graph_t *W)
{
    for (int64_t i = 0; i < W->n; i++) {
        for (int64_t k = W->colptr[i]; k < W->colptr[i + 1]; k++) {
            if (k > W->colptr[i] && W->row[k] <= W->row[k - 1]) {
                printf("row %ld not sorted\n", (long)i);
                return 1;
            }
            if (W->val[k] <= 0) {
                printf("nonpositive weight in row %ld\n", (long)i);
                return 1;
            }
            if (find_entry(W, W->row[k], i) < 0) {
                printf("not symmetric at (%ld, %ld)\n", (long)i, (long)W->row[k]);
                return 1;
            }
        }
    }
    return 0;
}
/* max_i |F_i - Fref_i| / max_i |Fref_i|, over both coordinates */
static double rel_err(const double *F, const double *Fref, int64_t n)
{
    double emax = 0, fmax = 0;
    for (int64_t k = 0; k < 2 * n; k++) {
        double e = fabs(F[k] - Fref[k]);
        if (e > emax)
            emax = e;
        if (fabs(Fref[k]) > fmax)
            fmax = fabs(Fref[k]);
    }
    return emax / fmax;
}
int main(int argc, char **argv)
{
       if (argc < 2) {
        printf("usage: %s data/<graph> [m]\n", argv[0]);
        return 1;
    }
    int m = (argc > 2) ? atoi(argv[2]) : 8;   /* leaf capacity */

    graph_t W;
    state_t S;
    ref_t   R;
    if (read_graph(argv[1], &W) != 0 || read_state(argv[1], &S) != 0 || read_ref(argv[1], &R) != 0) {
        printf("cannot read %s\n", argv[1]);
        return 1;
    }
    printf("%s: n = %ld, nnz = %ld, max degree = %ld\n",
           argv[1], (long)W.n, (long)W.nnz, (long)max_degree(&W));

    if (check_graph(&W) != 0)
        return 1;

    double *beta = malloc(W.n * sizeof(double));
    double *P    = malloc(W.nnz * sizeof(double));

    double t0 = now();
    compute_P_seq(&W, S.u, beta, P);
    double t1 = now();
    printf("P: %.3f s\n", t1 - t0);

    /* checks */
    double err = 0;
    for (int64_t i = 0; i < W.n; i++) {
        double e = fabs(beta[i] - R.beta[i]);
        if (e > err)
            err = e;
    }
    double sum = 0;
    for (int64_t k = 0; k < W.nnz; k++)
        sum += P[k];

    printf("max |beta - ref| = %.3e\n", err);
    printf("beta_1 = %.6f (ref %.6f)\n", beta[0], R.beta[0]);
    printf("sum P = %.12f (should be 1)\n", sum);
    
    /* Step 2: attraction */
    double *Fattr = malloc(2 * W.n * sizeof(double));

    t0 = now();
    compute_attr_seq(&W, P, S.Y0, Fattr);
    t1 = now();
    printf("attr: %.3f s\n", t1 - t0);
    printf("Fattr rel err = %.3e\n", rel_err(Fattr, R.Fattr, W.n));
    printf("Fattr_1 = (%.6e, %.6e) (ref (%.6e, %.6e))\n",
           Fattr[0], Fattr[1], R.Fattr[0], R.Fattr[1]);   

    /* Step 3: Barnes-Hut repulsion */
    double *Fbh = malloc(2 * W.n * sizeof(double));
    tree_t T;

    t0 = now();
    tree_build_seq(&T, W.n, S.Y0, m);
    t1 = now();
    tree_com_seq(&T, S.Y0);
    double t2 = now();
    double Zbh = compute_rep_bh_seq(&T, W.n, S.Y0, S.theta, Fbh);
    double t3 = now();

        printf("tree (m = %d): %ld cells, max depth %d\n", m, (long)T.ncells, T.maxdepth);
    printf("build: %.3f s, com: %.3f s, walk: %.3f s\n", t1 - t0, t2 - t1, t3 - t2);
        double t_bh = t3 - t0;
    printf("Z bh = %.2f (ref %.2f), rel err = %.3e\n", Zbh, R.Z, fabs(Zbh - R.Z) / R.Z);
    printf("Frep bh rel err = %.3e\n", rel_err(Fbh, R.Frep, W.n));

    /* one line for scripts: grep "^CSV_BH" */
    const char *name = strrchr(argv[1], '/');
    name = (name != NULL) ? name + 1 : argv[1];
    printf("CSV_BH,%s,%ld,%d,%ld,%d,%.6f,%.6f,%.6f,%.3e,%.3e\n",
           name, (long)W.n, m, (long)T.ncells, T.maxdepth,
           t1 - t0, t2 - t1, t3 - t2,
           fabs(Zbh - R.Z) / R.Z, rel_err(Fbh, R.Frep, W.n));

    /* Step 4: gradient */
    double *G    = malloc(2 * W.n * sizeof(double));
    double *Gref = malloc(2 * W.n * sizeof(double));

    t0 = now();
    compute_grad_seq(W.n, Fattr, Fbh, G);
    t1 = now();
    printf("grad: %.3f s\n", t1 - t0);

    compute_grad_seq(W.n, R.Fattr, R.Frep, Gref);
    printf("grad rel err = %.3e\n", rel_err(G, Gref, W.n));

    save_results(argv[1], "seq", W.n, Zbh, beta, Fattr, Fbh, G);
    /* Step 3 (check): direct O(n^2) repulsion, small graphs only */
    if (W.n <= 50000) {
        double *Fdir = malloc(2 * W.n * sizeof(double));

        t0 = now();
        double Zdir = compute_rep_direct_seq(W.n, S.Y0, Fdir);
        t1 = now();
        printf("direct rep: %.3f s\n", t1 - t0);
                double t_dir = t1 - t0;
        printf("Z direct = %.2f (ref %.2f), rel err = %.3e\n",
               Zdir, R.Z, fabs(Zdir - R.Z) / R.Z);
        printf("Frep direct rel err = %.3e\n", rel_err(Fdir, R.Frep, W.n));
        printf("Frep_1 = (%.6e, %.6e) (ref (%.6e, %.6e))\n",
               Fdir[0], Fdir[1], R.Frep[0], R.Frep[1]);
         printf("BH vs direct: Z rel err = %.3e, Frep rel err = %.3e\n",
               fabs(Zbh - Zdir) / Zdir, rel_err(Fbh, Fdir, W.n));
                       printf("BH speedup over direct: %.1fx\n", t_dir / t_bh);
        

        free(Fdir);
    } free(Fattr);
        free(Fbh);
            free(G);
    free(Gref);
    tree_free(&T);
    free(beta);
    free(P);
    free_graph(&W);
    return 0;
}