/*
 * Prints the number of blue edges in each t-subgraph for graph read from stdin
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#define NED_MAX 2048

void print_nbcount(int sp[], int nv, int t, int igraph)
{
    int c[t+2]; /* requires C99 auto dynamic arrays */
    int j, k, sum, nedc = t*(t-1)/2, isub = 0;

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

        printf("%8d %8d %8d\n", igraph, isub++, (sum + nedc)/2);

        /* FIND j */
        j = 0;
        while (c[j] + 1 == c[j+1]) { c[j] = j; j++; }

        /* DONE? */
        if (j == t) break;

        c[j]++;
    }
}

int main(int argc, char *argv[])
{
    int sp[NED_MAX];
    int nv, r, t, ned, i, igraph = 0;

    if (argc != 2)
    {
        fprintf(stderr, "Usage: %s t\n", argv[0]);
        fprintf(stderr, "Reads graph string from stdin\n");
        exit(EXIT_FAILURE);
    }

    t = atoi(argv[1]);

    printf("# %6s %8s %8s\n", "graph", "sub", "nb");

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

        print_nbcount(sp, nv, t, igraph++);
    }

    return EXIT_SUCCESS;
}
