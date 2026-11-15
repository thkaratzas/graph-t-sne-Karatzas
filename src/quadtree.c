#include <stdio.h>
#include <stdlib.h>
#include "tsne.h"

#define MAX_DEPTH 50   /* guard against many identical points */

/* 0: left-bottom, 1: right-bottom, 2: left-top, 3: right-top */
static int quadrant(const cell_t *c, double x, double y)
{
    return (x >= c->cx) + 2 * (y >= c->cy);
}

static cell_t *new_cell(double cx, double cy, double h, int64_t start, int64_t count)
{
    cell_t *c = malloc(sizeof(cell_t));
    c->cx = cx;
    c->cy = cy;
    c->h  = h;
    c->start = start;
    c->count = count;
    c->mx = 0;
    c->my = 0;
    c->leaf = 1;
    for (int q = 0; q < 4; q++)
        c->child[q] = NULL;
    return c;
}

/* ---------- shared pieces (all versions) ---------- */

/* Bounding square, identity perm, root cell.
   Returns a work buffer of size n; the caller frees it. */
int64_t *tree_init(tree_t *T, int64_t n, const double *Y, int m)
{
    double xmin = Y[0], xmax = Y[0], ymin = Y[1], ymax = Y[1];
    for (int64_t i = 1; i < n; i++) {
        double x = Y[2 * i], y = Y[2 * i + 1];
        if (x < xmin) xmin = x;
        if (x > xmax) xmax = x;
        if (y < ymin) ymin = y;
        if (y > ymax) ymax = y;
    }
    double cx = (xmin + xmax) / 2;
    double cy = (ymin + ymax) / 2;
    double h  = (xmax - xmin > ymax - ymin ? xmax - xmin : ymax - ymin) / 2;

    T->m = m;
    T->perm = malloc(n * sizeof(int64_t));
    for (int64_t i = 0; i < n; i++)
        T->perm[i] = i;

    T->root = new_cell(cx, cy, h, 0, n);
    return malloc(n * sizeof(int64_t));
}

/* Split cell c into its quadrants (one level, no recursion).
   Returns 0 if c stays a leaf, 1 if it got children. */
int tree_split(cell_t *c, int64_t *perm, int64_t *tmp, const double *Y, int m, int depth)
{
    if (c->count <= m || depth >= MAX_DEPTH)
        return 0;
    c->leaf = 0;

    int64_t *p = &perm[c->start];
    int64_t *t = &tmp[c->start];

    /* count points per quadrant */
    int64_t cnt[4] = {0, 0, 0, 0};
    for (int64_t k = 0; k < c->count; k++)
        cnt[quadrant(c, Y[2 * p[k]], Y[2 * p[k] + 1])]++;

    /* where each quadrant starts */
    int64_t off[4];
    off[0] = 0;
    for (int q = 1; q < 4; q++)
        off[q] = off[q - 1] + cnt[q - 1];

    /* group the points by quadrant: into tmp, then back */
    int64_t pos[4] = {off[0], off[1], off[2], off[3]};
    for (int64_t k = 0; k < c->count; k++) {
        int q = quadrant(c, Y[2 * p[k]], Y[2 * p[k] + 1]);
        t[pos[q]++] = p[k];
    }
    for (int64_t k = 0; k < c->count; k++)
        p[k] = t[k];

    /* children */
    double h2 = c->h / 2;
    for (int q = 0; q < 4; q++) {
        if (cnt[q] == 0)
            continue;
        double cx = c->cx + ((q & 1) ? h2 : -h2);
        double cy = c->cy + ((q & 2) ? h2 : -h2);
        c->child[q] = new_cell(cx, cy, h2, c->start + off[q], cnt[q]);
    }
    return 1;
}

static int64_t count_cells(const cell_t *c)
{
    int64_t k = 1;
    for (int q = 0; q < 4; q++)
        if (c->child[q])
            k += count_cells(c->child[q]);
    return k;
}

/* depth of the deepest leaf below c (c itself at depth 0) */
static int max_depth(const cell_t *c)
{
    int d = 0;
    for (int q = 0; q < 4; q++) {
        if (c->child[q]) {
            int dq = 1 + max_depth(c->child[q]);
            if (dq > d)
                d = dq;
        }
    }
    return d;
}

void tree_finish(tree_t *T)
{
    T->ncells = count_cells(T->root);
    T->maxdepth = max_depth(T->root);
}

/* Center of mass of one cell: from its points (leaf)
   or from its children, which must be done already. */
void tree_com_cell(cell_t *c, const int64_t *perm, const double *Y)
{
    double sx = 0, sy = 0;
    if (c->leaf) {
        for (int64_t k = c->start; k < c->start + c->count; k++) {
            sx += Y[2 * perm[k]];
            sy += Y[2 * perm[k] + 1];
        }
    } else {
        for (int q = 0; q < 4; q++) {
            cell_t *ch = c->child[q];
            if (ch == NULL)
                continue;
            sx += ch->count * ch->mx;
            sy += ch->count * ch->my;
        }
    }
    c->mx = sx / c->count;
    c->my = sy / c->count;
}

/* ---------- sequential version ---------- */

static void build(cell_t *c, int64_t *perm, int64_t *tmp, const double *Y, int m, int depth)
{
    if (!tree_split(c, perm, tmp, Y, m, depth))
        return;
    for (int q = 0; q < 4; q++)
        if (c->child[q])
            build(c->child[q], perm, tmp, Y, m, depth + 1);
}

void tree_build_seq(tree_t *T, int64_t n, const double *Y, int m)
{
    int64_t *tmp = tree_init(T, n, Y, m);
    build(T->root, T->perm, tmp, Y, m, 0);
    free(tmp);
    tree_finish(T);
}

static void com(cell_t *c, const int64_t *perm, const double *Y)
{
    if (!c->leaf)
        for (int q = 0; q < 4; q++)
            if (c->child[q])
                com(c->child[q], perm, Y);
    tree_com_cell(c, perm, Y);
}

void tree_com_seq(tree_t *T, const double *Y)
{
    com(T->root, T->perm, Y);
}

/* ---------- diagnostics and cleanup ---------- */

/* Leaves at the depth limit: they hold points that could not be separated. */
static void deep_leaves(const cell_t *c, int depth,
                        int64_t *nleaves, int64_t *npoints, int64_t *maxcount)
{
    if (c->leaf) {
        if (depth >= MAX_DEPTH) {
            (*nleaves)++;
            *npoints += c->count;
            if (c->count > *maxcount)
                *maxcount = c->count;
        }
        return;
    }
    for (int q = 0; q < 4; q++)
        if (c->child[q])
            deep_leaves(c->child[q], depth + 1, nleaves, npoints, maxcount);
}

void tree_deep_report(const tree_t *T)
{
    int64_t nleaves = 0, npoints = 0, maxcount = 0;
    deep_leaves(T->root, 0, &nleaves, &npoints, &maxcount);
    printf("leaves at depth %d: %ld, points in them: %ld, largest: %ld\n",
           MAX_DEPTH, (long)nleaves, (long)npoints, (long)maxcount);
}

static void free_cell(cell_t *c)
{
    for (int q = 0; q < 4; q++)
        if (c->child[q])
            free_cell(c->child[q]);
    free(c);
}

void tree_free(tree_t *T)
{
    free_cell(T->root);
    free(T->perm);
}