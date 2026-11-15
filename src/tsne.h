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
    int      maxdepth;      /* depth of the deepest leaf (root = 0) */
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
void tree_deep_report(const tree_t *T);

/* shared tree pieces (quadtree.c) */
int64_t *tree_init(tree_t *T, int64_t n, const double *Y, int m);
int  tree_split(cell_t *c, int64_t *perm, int64_t *tmp, const double *Y, int m, int depth);
void tree_finish(tree_t *T);
void tree_com_cell(cell_t *c, const int64_t *perm, const double *Y);
double compute_rep_bh_seq(const tree_t *T, int64_t n, const double *Y, double theta, double *F);

/* Step 4: gradient G[2n] = 4 (Fattr - Frep) */
void compute_grad_seq(int64_t n, const double *Fattr, const double *Frep, double *G);

/* results/<graph>-<version>.bin (io.c) */
int save_results(const char *base, const char *version, int64_t n, double Z,
                 const double *beta, const double *Fattr, const double *Frep, const double *G);


/* OpenMP versions (*_omp.c), same arguments as the sequential ones */
void   compute_P_omp(const graph_t *W, double u, double *beta, double *P);
void   compute_attr_omp(const graph_t *W, const double *P, const double *Y, double *F);
double compute_rep_direct_omp(int64_t n, const double *Y, double *F);
void   tree_build_omp(tree_t *T, int64_t n, const double *Y, int m);
void   tree_com_omp(tree_t *T, const double *Y);
double compute_rep_bh_omp(const tree_t *T, int64_t n, const double *Y, double theta, double *F);
void   compute_grad_omp(int64_t n, const double *Fattr, const double *Frep, double *G);
/* helpers (p_common.c) */
double row_beta(const double *w, int64_t d, double u, double *p, double *lr);
int64_t find_entry(const graph_t *W, int64_t j, int64_t i);
int64_t max_degree(const graph_t *W);

/* Barnes-Hut walk for one point (bh_common.c) */
void bh_walk(const cell_t *c, const int64_t *perm, const double *Y,
             int64_t i, double theta, double *acc);

#endif