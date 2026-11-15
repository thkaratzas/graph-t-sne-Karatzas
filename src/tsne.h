#ifndef TSNE_H
#define TSNE_H
#include <stdint.h>
#include "data.h"

/* Step 1: beta[n] and P[nnz] (P has the same sparsity as W) */
void compute_P_seq(const graph_t *W, double u, double *beta, double *P);

/* Step 2: attraction F[2n] at points Y[2n] */
void compute_attr_seq(const graph_t *W, const double *P, const double *Y, double *F);

/* helpers (p_common.c) */
double row_beta(const double *w, int64_t d, double u, double *p, double *lr);
int64_t find_entry(const graph_t *W, int64_t j, int64_t i);
int64_t max_degree(const graph_t *W);

#endif