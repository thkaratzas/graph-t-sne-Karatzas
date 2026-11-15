#include <stdlib.h>
#include <omp.h>
#include "tsne.h"

/* A child with fewer points than this is done by the same thread,
   without a new task: a task costs more than a small subtree. */
#define TASK_CUTOFF 10000

static void build_omp(cell_t *c, int64_t *perm, int64_t *tmp, const double *Y, int m, int depth)
{
    if (!tree_split(c, perm, tmp, Y, m, depth))
        return;

    for (int q = 0; q < 4; q++) {
        cell_t *ch = c->child[q];
        if (ch == NULL)
            continue;
        if (ch->count > TASK_CUTOFF) {
            #pragma omp task
            build_omp(ch, perm, tmp, Y, m, depth + 1);
        } else {
            build_omp(ch, perm, tmp, Y, m, depth + 1);
        }
    }
    #pragma omp taskwait
}

void tree_build_omp(tree_t *T, int64_t n, const double *Y, int m)
{
    int64_t *tmp = tree_init(T, n, Y, m);

    #pragma omp parallel
    #pragma omp single
    build_omp(T->root, T->perm, tmp, Y, m, 0);

    free(tmp);
    tree_finish(T);
}

static void com_omp(cell_t *c, const int64_t *perm, const double *Y)
{
    if (!c->leaf) {
        for (int q = 0; q < 4; q++) {
            cell_t *ch = c->child[q];
            if (ch == NULL)
                continue;
            if (ch->count > TASK_CUTOFF) {
                #pragma omp task
                com_omp(ch, perm, Y);
            } else {
                com_omp(ch, perm, Y);
            }
        }
        #pragma omp taskwait   /* children first, then this cell */
    }
    tree_com_cell(c, perm, Y);
}

void tree_com_omp(tree_t *T, const double *Y)
{
    #pragma omp parallel
    #pragma omp single
    com_omp(T->root, T->perm, Y);
}