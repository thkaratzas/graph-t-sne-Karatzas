#include <stdio.h>
#include <string.h>
#include "tsne.h"

/* File layout (all little-endian, as in memory):
     int64  n
     double Z
     double beta[n]
     double Fattr[2n]   point i at [2i], [2i+1]
     double Frep[2n]
     double G[2n]       gradient 4 (Fattr - Frep)                          */
int save_results(const char *base, const char *version, int64_t n, double Z,
                 const double *beta, const double *Fattr, const double *Frep, const double *G)
{
    /* graph name: the part of base after the last '/' */
    const char *name = strrchr(base, '/');
    name = (name != NULL) ? name + 1 : base;

    char path[512];
    snprintf(path, sizeof(path), "results/%s-%s.bin", name, version);

    FILE *f = fopen(path, "wb");
    if (f == NULL) {
        printf("cannot write %s\n", path);
        return 1;
    }
    fwrite(&n, sizeof(int64_t), 1, f);
    fwrite(&Z, sizeof(double), 1, f);
    fwrite(beta,  sizeof(double), n, f);
    fwrite(Fattr, sizeof(double), 2 * n, f);
    fwrite(Frep,  sizeof(double), 2 * n, f);
    fwrite(G,     sizeof(double), 2 * n, f);
    fclose(f);

    printf("saved %s\n", path);
    return 0;
}