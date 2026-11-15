#include <stdio.h>
#include <pthread.h>
#include <omp.h>

void *pt_hello(void *arg)
{
    int id = *(int *)arg;
    printf("pthread %d\n", id);
    return NULL;
}

int main(void)
{
    pthread_t t[2];
    int ids[2] = {0, 1};
    for (int i = 0; i < 2; i++)
        pthread_create(&t[i], NULL, pt_hello, &ids[i]);
    for (int i = 0; i < 2; i++)
        pthread_join(t[i], NULL);

    #pragma omp parallel
    printf("omp thread %d of %d\n", omp_get_thread_num(), omp_get_num_threads());

    return 0;
}
