#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#define NED_MAX 2048

int clique_count(int sp[], int nv, int color, int n)
{
    int c[n+2]; /* uses C99 auto dynamic arrays */
    int nedc, count, j, k, sum;

    nedc = n*(n-1)/2;
    count = 0;

    /* INITIALIZE */
    c[n] = nv;
    c[n+1] = 0;
    for (j = 0; j < n; j++) c[j] = j;

    while (1)
    {
        sum = 0;

        /* iterate over edges in this subgraph */
        for (k = 0; k < n; k++)
            for (j = 0; j < k; j++)
                sum += sp[c[k]*(c[k]-1)/2 + c[j]];

        if (sum == color*nedc) count++;

        /* FIND j */
        j = 0;
        while (c[j] + 1 == c[j+1]) { c[j] = j; j++; }

        /* DONE? */
        if (j == n) break;

        c[j]++;
    }

    return count;
}

int energy(int sp[], int nv, int r, int s)
{
    return clique_count(sp, nv, -1, r) + clique_count(sp, nv, 1, s);
}

int main(int argc, char *argv[])
{
    int sp[NED_MAX];
    int nv, r, s, ned, i;

    if (argc != 3)
    {
        fprintf(stderr, "Usage: %s r s\n", argv[0]);
        fprintf(stderr, "Reads graph string from stdin\n");
        exit(EXIT_FAILURE);
    }

    r = atoi(argv[1]);
    s = atoi(argv[2]);

    while (scanf("%d", &nv) > 0)
    {
        ned = nv*(nv-1)/2;
        assert(ned < NED_MAX);

        for (i = 0; i < ned; i++)
        {
            if (scanf("%d", &sp[i]) <= 0)
            {
                fprintf(stderr, "error in input\n");
                return EXIT_FAILURE;
            }

            if (sp[i] == 0) sp[i] = -1;
        }

        printf("%d\n", energy(sp, nv, r, s));
    }

    return EXIT_SUCCESS;
}
