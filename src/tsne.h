#ifndef TSNE_H
#define TSNE_H
#include <stdint.h>
#include "data.h"

/* Quadtree cell. Its points are perm[start .. start+count-1]. */
typedef struct cell {
    double cx, cy;          /* center of the square */
    double h;               /* half side */
    int64_t start, count;   /* points of this cell in perm[] */
    double mx, my;          /* center of mass */
    struct cell *child[4];  /* NULL if empty or leaf */
    int leaf;
} cell_t;

typedef struct {
    cell_t  *root;
    int64_t *perm;          /* point indices, grouped by cell */
    int      m;             /* leaf capacity */
    int64_t  ncells;
} tree_t;
/* Step 1: beta[n] and P[nnz] (P has the same sparsity as W) */
void compute_P_seq(const graph_t *W, double u, double *beta, double *P);

/* Step 2: attraction F[2n] at points Y[2n] */
void compute_attr_seq(const graph_t *W, const double *P, const double *Y, double *F);

/* Step 3 (check): direct O(n^2) repulsion F[2n], returns Z */
double compute_rep_direct_seq(int64_t n, const double *Y, double *F);

/* Step 3: quadtree (quadtree.c) and Barnes-Hut (bh_seq.c) */
void tree_build_seq(tree_t *T, int64_t n, const double *Y, int m);
void tree_com_seq(tree_t *T, const double *Y);
void tree_free(tree_t *T);
double compute_rep_bh_seq(const tree_t *T, int64_t n, const double *Y, double theta, double *F);
/* helpers (p_common.c) */
double row_beta(const double *w, int64_t d, double u, double *p, double *lr);
int64_t find_entry(const graph_t *W, int64_t j, int64_t i);
int64_t max_degree(const graph_t *W);

#endif