#include <stdio.h>
#include <cilk/cilk.h>

int main(void)
{
    long sum = 0;
    int a[1000];
    cilk_for (int i = 0; i < 1000; i++)
        a[i] = i;
    for (int i = 0; i < 1000; i++)
        sum += a[i];
    printf("cilk sum = %ld (should be 499500)\n", sum);
    return 0;
}
