#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>
#include <stdatomic.h>
#include "tsne.h"

/* settings, read once by pt_init() */
static int     pt_p       = 1;     /* number of threads */
static int     pt_dynamic = 1;     /* 1: dynamic, 0: static */
static int64_t pt_chunk   = 256;   /* chunk size for dynamic */

void pt_init(void)
{
    char *s = getenv("PT_NUM_THREADS");
    pt_p = (s != NULL) ? atoi(s) : (int)sysconf(_SC_NPROCESSORS_ONLN);
    if (pt_p < 1)
        pt_p = 1;

    s = getenv("PT_SCHEDULE");
    pt_dynamic = !(s != NULL && strcmp(s, "static") == 0);

    s = getenv("PT_CHUNK");
    if (s != NULL)
        pt_chunk = atoll(s);
    if (pt_chunk < 1)
        pt_chunk = 1;
}

int pt_nthreads(void)            { return pt_p; }
const char *pt_schedule(void)    { return pt_dynamic ? "dynamic" : "static"; }
int64_t pt_get_chunk(void)       { return pt_chunk; }

/* what all threads of one pt_for share */
typedef struct {
    int64_t n;
    pt_body_t body;
    void *arg;
    atomic_llong next;   /* dynamic: first iteration not yet taken */
} shared_t;

/* what each thread gets */
typedef struct {
    shared_t *s;
    int tid;
} worker_t;

static void *worker(void *v)
{
    worker_t *w = v;
    shared_t *s = w->s;

    if (pt_dynamic) {
        /* take the next chunk until none is left */
        int64_t lo;
        while ((lo = atomic_fetch_add(&s->next, pt_chunk)) < s->n) {
            int64_t hi = (lo + pt_chunk < s->n) ? lo + pt_chunk : s->n;
            s->body(lo, hi, w->tid, s->arg);
        }
    } else {
        /* one contiguous block per thread */
        int64_t lo = s->n * w->tid / pt_p;
        int64_t hi = s->n * (w->tid + 1) / pt_p;
        if (lo < hi)
            s->body(lo, hi, w->tid, s->arg);
    }
    return NULL;
}

/* Run body(lo, hi, tid, arg) over [0, n), split among pt_p threads.
   The calling thread works too, as thread 0. Returns when all is done. */
void pt_for(int64_t n, pt_body_t body, void *arg)
{
    shared_t s;
    s.n = n;
    s.body = body;
    s.arg = arg;
    atomic_init(&s.next, 0);

    pthread_t th[pt_p];
    worker_t  w[pt_p];
    for (int t = 0; t < pt_p; t++) {
        w[t].s = &s;
        w[t].tid = t;
    }

    for (int t = 1; t < pt_p; t++)
        pthread_create(&th[t], NULL, worker, &w[t]);
    worker(&w[0]);
    for (int t = 1; t < pt_p; t++)
        pthread_join(th[t], NULL);
}