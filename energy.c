#include <stdio.h>
#include <stdlib.h>

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
    FILE *fp;
    int *sp;
    int nv, r, s, ned, i = 0;

    if (argc != 4)
    {
        fprintf(stderr, "Usage: %s graph_file r s\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    fp = fopen(argv[1], "r");
    r = atoi(argv[2]);
    s = atoi(argv[3]);

    fscanf(fp, "%d", &nv);
    ned = nv*(nv-1)/2;
    sp = (int*) malloc(ned * sizeof(int));

    while (fscanf(fp, "%d", &sp[i]) && i < ned)
    {
        if (sp[i] == 0) sp[i] = -1;
        i++;
    }

    printf("%d\n", energy(sp, nv, r, s));

    return EXIT_SUCCESS;
}
