/* Compare two result files:  ./bin/compare results/A.bin results/B.bin
   Prints max_k |a_k - b_k| / max_k |b_k| for each array, and |Za - Zb| / |Zb|. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

typedef struct {
    int64_t n;
    double Z;
    double *beta, *Fattr, *Frep, *G;
} res_t;

static int load(const char *path, res_t *r)
{
    FILE *f = fopen(path, "rb");
    if (f == NULL) {
        printf("cannot open %s\n", path);
        return 1;
    }
    if (fread(&r->n, sizeof(int64_t), 1, f) != 1 || fread(&r->Z, sizeof(double), 1, f) != 1) {
        printf("cannot read header of %s\n", path);
        fclose(f);
        return 1;
    }
    int64_t n = r->n;
    r->beta  = malloc(n * sizeof(double));
    r->Fattr = malloc(2 * n * sizeof(double));
    r->Frep  = malloc(2 * n * sizeof(double));
    r->G     = malloc(2 * n * sizeof(double));
    size_t got = fread(r->beta, sizeof(double), n, f)
               + fread(r->Fattr, sizeof(double), 2 * n, f)
               + fread(r->Frep, sizeof(double), 2 * n, f)
               + fread(r->G, sizeof(double), 2 * n, f);
    fclose(f);
    if (got != (size_t)(7 * n)) {
        printf("%s is too short\n", path);
        return 1;
    }
    return 0;
}

static double rel(const double *a, const double *b, int64_t len)
{
    double emax = 0, bmax = 0;
    for (int64_t k = 0; k < len; k++) {
        double e = fabs(a[k] - b[k]);
        if (e > emax)
            emax = e;
        if (fabs(b[k]) > bmax)
            bmax = fabs(b[k]);
    }
    return (bmax > 0) ? emax / bmax : emax;
}

int main(int argc, char **argv)
{
    if (argc < 3) {
        printf("usage: %s A.bin B.bin\n", argv[0]);
        return 1;
    }
    res_t a, b;
    if (load(argv[1], &a) != 0 || load(argv[2], &b) != 0)
        return 1;
    if (a.n != b.n) {
        printf("different n: %ld vs %ld\n", (long)a.n, (long)b.n);
        return 1;
    }

    int64_t n = a.n;
    printf("%s vs %s (n = %ld)\n", argv[1], argv[2], (long)n);
    printf("  Z     %.3e\n", fabs(a.Z - b.Z) / fabs(b.Z));
    printf("  beta  %.3e\n", rel(a.beta, b.beta, n));
    printf("  Fattr %.3e\n", rel(a.Fattr, b.Fattr, 2 * n));
    printf("  Frep  %.3e\n", rel(a.Frep, b.Frep, 2 * n));
    printf("  grad  %.3e\n", rel(a.G, b.G, 2 * n));
    return 0;
}