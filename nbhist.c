/*
 * Prints a histogram of the number of blue edges in each t-subgraph for graph
 * read from stdin
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NBIN_MAX    50
#define NED_MAX     2000
#define NEDC_MAX    50      /* max number of edges in a clique (restricts max
                               value of argument t) */

void accumulate(double count_sum[], double count2_sum[], int sp[], int nv, int t, int nedc)
{
    int count[NBIN_MAX];
    int c[NEDC_MAX+2];
    int j, k, sum, nb;

    for (nb = 0; nb < NBIN_MAX; nb++) count[nb] = 0;

    /* INITIALIZE */
    c[t] = nv;
    c[t+1] = 0;
    for (j = 0; j < t; j++) c[j] = j;

    for (;;)
    {
        sum = 0;

        /* iterate over edges in this subgraph */
        for (k = 0; k < t; k++)
            for (j = 0; j < k; j++)
                sum += sp[c[k]*(c[k]-1)/2 + c[j]];

        nb = (sum + nedc)/2;
#ifdef DEBUG
        assert(nb < NBIN_MAX);
#endif
        count[nb]++;

        /* FIND j */
        j = 0;
        while (c[j] + 1 == c[j+1]) { c[j] = j; j++; }

        /* DONE? */
        if (j == t) break;

        c[j]++;
    }

    for (nb = 0; nb < NBIN_MAX; nb++)
    {
        count_sum[nb] += (double) count[nb];
        count2_sum[nb] += (double) count[nb]*count[nb];
    }
}

int main(int argc, char *argv[])
{
    double count_sum[NBIN_MAX], count2_sum[NBIN_MAX];
    int sp[NED_MAX];
    double var;
    int i, nv, ned, t, nedc, ngraph = 0, imax = NBIN_MAX-1;

    if (argc != 2)
    {
        fprintf(stderr, "Usage: %s t\n", argv[0]);
        fprintf(stderr, "Reads graph string from stdin\n");
        exit(EXIT_FAILURE);
    }

    t = atoi(argv[1]);
    nedc = t*(t-1)/2;
    assert(nedc <= NEDC_MAX);

    for (i = 0; i < NBIN_MAX; i++) { count_sum[i] = 0.; count2_sum[i] = 0.; }

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

        accumulate(count_sum, count2_sum, sp, nv, t, nedc);
        ngraph += 1;
    }

    while (count_sum[imax] == 0.) imax--;
    i = 0; while (count_sum[i] == 0.) i++;

    if (ngraph > 0)
    {
        if (ngraph > 1)
        {
            printf("# %2s %10s %10s\n", "nb", "count_av", "err");
            for (; i <= imax; i++)
            {
                var = (count2_sum[i] - pow(count_sum[i],2)/ngraph) / ngraph;
                printf("%4d %10d %10d\n",
                    i, (int) count_sum[i]/ngraph, (int) sqrt(var/(ngraph-1.)));
            }
        }
        else    /* only one graph in the input */
        {
            printf("# %2s %10s\n", "nb", "count");
            for (; i <= imax; i++)
                printf("%4d %10d\n", i, (int) count_sum[i]);
        }
    }

    return EXIT_SUCCESS;
}
