/* Does a brute-force computation of the "soft" energy used in ramsey2.c.
 * Useful for making sure updates are working properly */

double part_energy(int sp[NED], int color, int t, double et[])
{
    int c[S+2];
    int count, j, k;
    double en = 0.;

    /* INITIALIZE */
    c[t] = NV;
    c[t+1] = 0;
    for (j = 0; j < t; j++) c[j] = j;

    while (1)
    {
        count = 0;

        /* iterate over edges in this subgraph */
        for (k = 0; k < t; k++)
            for (j = 0; j < k; j++)
                if (sp[c[k]*(c[k]-1)/2 + c[j]] == color)
                    count++;

        en += et[count];

        /* FIND j */
        j = 0;
        while (c[j] + 1 == c[j+1]) { c[j] = j; j++; }

        /* DONE? */
        if (j == t) break;

        c[j]++;
    }

    return en;
}

double debug_energy(int sp[NED], double er[NEDR+1], double es[NEDS+1])
{
    return part_energy(sp, 1, R, er) + part_energy(sp, -1, S, es);
}
