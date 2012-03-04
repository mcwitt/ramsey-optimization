#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ramsey.h"
#include "sga.h"

void decode(SGA_allele_t chrom[], int sp[])
{
    int i;

    for (i = 0; i < NED; i++)
        sp[i] = (chrom[i] == 1) ? 1 : -1;
}

double objfunc(SGA_allele_t chrom[])
{
    int sp[NED];

    decode(chrom, sp);
    return -R_energy(sp);
}

int main(int argc, char *argv[])
{
    SGA_t sga;
    char filename[256];
    double pcross, pmutate;
    int popsize;
    int sp[NED];
    int igen, ngen, seed;
    int ncross, nmutation;

    if (argc != 6)
    {
        fprintf(stderr, "Usage: %s ngen popsize pcross pmutate seed\n", argv[0]);
        fprintf(stderr, "Compiled for (%d, %d, %d)\n", R, S, NV);
        exit(EXIT_FAILURE);
    }

    ngen    = atoi(argv[1]);
    popsize = atoi(argv[2]);
    pcross  = atof(argv[3]);
    pmutate = atof(argv[4]);
    seed    = atoi(argv[5]);

    R_init(seed);
    SGA_init(&sga, popsize, NED, objfunc, pcross, pmutate, seed);

    sprintf(filename, "%d-%d-%d_%d.graph", R, S, NV, seed);
    printf("#%8s %9s %9s %9s %9s %9s %9s %9s\n",
            "gen", "emin", "favg", "fvar", "fmin", "fmax", "ncross", "nmutation");

    for (igen = 0; igen < ngen; igen++)
    {
        SGA_advance(&sga, &ncross, &nmutation);

        if (igen % 10 == 0)
        {
            printf("%9d %9.3g %9.3g %9.3g %9.3g %9.3g %9d %9d\n",
                    igen,
                    -sga.objective[sga.fittest],
                    sga.favg,
                    sga.fvar,
                    sga.fmin,
                    sga.fmax,
                    ncross,
                    nmutation
                  );
            fflush(stdout);
        }

        if (-sga.objective[sga.fittest] < 10e-9) break;
    }

    decode(sga.chrom[sga.fittest], sp);
    R_save_graph(sp, filename);

    return EXIT_SUCCESS;
}
