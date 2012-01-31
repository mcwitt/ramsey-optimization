#ifdef R
#if (R > S)
#define M R
#else
#define M S
#endif
#else
#define R S
#endif

int debug_energy(int sp[NED])
{
    int c[R+S+2];
    int j, k;
    int edgesum, energy;
    int nedr, neds;

    nedr = R*(R-1)/2;
    neds = S*(S-1)/2;
    energy = 0;

    /* INITIALIZE */
    c[R] = NV;
    c[R+1] = 0;
    for (j = 0; j < R; j++) c[j] = j;

    while (1)
    {
        edgesum = 0;

        /* iterate over edges in this subgraph */
        for (k = 0; k < R; k++)
            for (j = 0; j < k; j++)
                edgesum += sp[c[k]*(c[k]-1)/2 + c[j]];

        if (edgesum == -nedr) energy++;

        /* FIND j */
        j = 0;
        while (c[j] + 1 == c[j+1]) { c[j] = j; j++; }

        /* DONE? */
        if (j == R) break;

        c[j]++;
    }

    /* INITIALIZE */
    c[S] = NV;
    c[S+1] = 0;
    for (j = 0; j < S; j++) c[j] = j;

    while (1)
    {
        edgesum = 0;

        /* iterate over edges in this subgraph */
        for (k = 0; k < S; k++)
            for (j = 0; j < k; j++)
                edgesum += sp[c[k]*(c[k]-1)/2 + c[j]];

        if (edgesum == neds) energy++;

        /* FIND j */
        j = 0;
        while (c[j] + 1 == c[j+1]) { c[j] = j; j++; }

        /* DONE? */
        if (j == S) break;

        c[j]++;
    }

    return energy;
}
