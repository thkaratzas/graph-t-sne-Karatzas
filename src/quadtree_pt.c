#include <stdlib.h>
#include <pthread.h>
#include "tsne.h"

#define SPAWN_CUTOFF 10000   /* same role as TASK_CUTOFF / SPAWN_CUTOFF */

/* smallest D with 4^D >= p */
static int spawn_depth(int p)
{
    int d = 0, k = 1;
    while (k < p) {
        k *= 4;
        d++;
    }
    return d;
}

/* ---------- tree build ---------- */

typedef struct {
    cell_t *c;
    int64_t *perm, *tmp;
    const double *Y;
    int m, depth, maxspawn;
} build_args_t;

static void *build_pt(void *v)
{
    build_args_t *a = v;
    if (!tree_split(a->c, a->perm, a->tmp, a->Y, a->m, a->depth))
        return NULL;

    build_args_t ca[4];
    pthread_t th[4];
    int started[4] = {0, 0, 0, 0};

    for (int q = 0; q < 4; q++) {
        cell_t *ch = a->c->child[q];
        if (ch == NULL)
            continue;
        ca[q] = *a;              /* same arrays and settings */
        ca[q].c = ch;
        ca[q].depth = a->depth + 1;
        if (a->depth < a->maxspawn && ch->count > SPAWN_CUTOFF) {
            pthread_create(&th[q], NULL, build_pt, &ca[q]);
            started[q] = 1;
        } else {
            build_pt(&ca[q]);
        }
    }
    for (int q = 0; q < 4; q++)
        if (started[q])
            pthread_join(th[q], NULL);
    return NULL;
}

void tree_build_pt(tree_t *T, int64_t n, const double *Y, int m)
{
    int64_t *tmp = tree_init(T, n, Y, m);

    build_args_t a;
    a.c = T->root;
    a.perm = T->perm;
    a.tmp = tmp;
    a.Y = Y;
    a.m = m;
    a.depth = 0;
    a.maxspawn = spawn_depth(pt_nthreads());
    build_pt(&a);

    free(tmp);
    tree_finish(T);
}

/* ---------- centers of mass ---------- */

typedef struct {
    cell_t *c;
    const int64_t *perm;
    const double *Y;
    int depth, maxspawn;
} com_args_t;

static void *com_pt(void *v)
{
    com_args_t *a = v;
    cell_t *c = a->c;

    if (!c->leaf) {
        com_args_t ca[4];
        pthread_t th[4];
        int started[4] = {0, 0, 0, 0};

        for (int q = 0; q < 4; q++) {
            cell_t *ch = c->child[q];
            if (ch == NULL)
                continue;
            ca[q] = *a;
            ca[q].c = ch;
            ca[q].depth = a->depth + 1;
            if (a->depth < a->maxspawn && ch->count > SPAWN_CUTOFF) {
                pthread_create(&th[q], NULL, com_pt, &ca[q]);
                started[q] = 1;
            } else {
                com_pt(&ca[q]);
            }
        }
        for (int q = 0; q < 4; q++)   /* children first, then this cell */
            if (started[q])
                pthread_join(th[q], NULL);
    }
    tree_com_cell(c, a->perm, a->Y);
    return NULL;
}

void tree_com_pt(tree_t *T, const double *Y)
{
    com_args_t a;
    a.c = T->root;
    a.perm = T->perm;
    a.Y = Y;
    a.depth = 0;
    a.maxspawn = spawn_depth(pt_nthreads());
    com_pt(&a);
}