#include <stdlib.h>
#include <cilk/cilk.h>
#include "tsne.h"

/* A child with fewer points than this is done in the same strand,
   without a spawn: a spawn costs more than a small subtree. */
#define SPAWN_CUTOFF 10000

static void build_cilk(cell_t *c, int64_t *perm, int64_t *tmp, const double *Y, int m, int depth)
{
    if (!tree_split(c, perm, tmp, Y, m, depth))
        return;

    for (int q = 0; q < 4; q++) {
        cell_t *ch = c->child[q];
        if (ch == NULL)
            continue;
        if (ch->count > SPAWN_CUTOFF)
            cilk_spawn build_cilk(ch, perm, tmp, Y, m, depth + 1);
        else
            build_cilk(ch, perm, tmp, Y, m, depth + 1);
    }
    cilk_sync;
}

void tree_build_cilk(tree_t *T, int64_t n, const double *Y, int m)
{
    int64_t *tmp = tree_init(T, n, Y, m);
    build_cilk(T->root, T->perm, tmp, Y, m, 0);
    free(tmp);
    tree_finish(T);
}

static void com_cilk(cell_t *c, const int64_t *perm, const double *Y)
{
    if (!c->leaf) {
        for (int q = 0; q < 4; q++) {
            cell_t *ch = c->child[q];
            if (ch == NULL)
                continue;
            if (ch->count > SPAWN_CUTOFF)
                cilk_spawn com_cilk(ch, perm, Y);
            else
                com_cilk(ch, perm, Y);
        }
        cilk_sync;   /* children first, then this cell */
    }
    tree_com_cell(c, perm, Y);
}

void tree_com_cilk(tree_t *T, const double *Y)
{
    com_cilk(T->root, T->perm, Y);
}