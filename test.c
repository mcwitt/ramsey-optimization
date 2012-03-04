#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "ramsey.h"

int clique_count(int sp[], int nv, int color, int t)
{
    int c[t+2]; /* requires C99 auto dynamic arrays */
    int sum, j, k, sumc = color * t*(t-1)/2, count = 0;

    /* INITIALIZE */
    c[t] = nv;
    c[t+1] = 0;
    for (j = 0; j < t; j++) c[j] = j;

    while (1)
    {
        sum = 0;

        /* iterate over edges in this subgraph */
        for (k = 0; k < t; k++)
            for (j = 0; j < k; j++)
                sum += sp[c[k]*(c[k]-1)/2 + c[j]];

        if (sum == sumc) count++;

        /* FIND j */
        j = 0;
        while (c[j] + 1 == c[j+1]) { c[j] = j; j++; }

        /* DONE? */
        if (j == t) break;

        c[j]++;
    }

    return count;
}

int energy(int sp[], int nv, int r, int s)
{
    return clique_count(sp, nv, -1, r) + clique_count(sp, nv, 1, s);
}

int main()
{
    R_replica_t r;
    int i;
    uint32_t seed = 0;

    R_init(seed);

    for (i = 0; i < 50; i++)
    {
        R_init_replica(&r);
        R_randomize(&r, 0.5, 0);
        //printf("%d\n", energy(r.sp, NV, R, S));
        printf("%d\n", R_energy(r.sp));
    }

    return EXIT_SUCCESS;
}
