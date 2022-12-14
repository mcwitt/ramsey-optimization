#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ramsey.h"
#include "ga_r.h"

void decode(GA_allele_t chrom[], int sp[])
{
    int i;

    for (i = 0; i < NED; i++)
        sp[i] = (chrom[i] == 1) ? 1 : -1;
}

//double fitfunc(double energy) { return 1./(energy + 1.); }
double fitfunc(double energy) { return -energy; }

double objfunc(GA_allele_t chrom[])
{
    int sp[NED];

    decode(chrom, sp);
    return R_energy(sp);
}

int main(int argc, char *argv[])
{
    GA_t ga;
    char filename[64];
    double pcross, pmutate;
    int sp[NED];
    int popsize, igen, ngen, ncross, nmutation, seed;

    if (argc != 6)
    {
        fprintf(stderr, "Usage: %s popsize pcross pmutate ngen seed\n", argv[0]);
        fprintf(stderr, "Compiled for (%d, %d, %d)\n", R, S, NV);
        exit(EXIT_FAILURE);
    }

    popsize = atoi(argv[1]);
    pcross  = atof(argv[2]);
    pmutate = atof(argv[3]);
    ngen    = atoi(argv[4]);
    seed    = atoi(argv[5]);

    R_init(seed);
    GA_init(&ga, popsize, NED, objfunc, fitfunc, pcross, pmutate, seed);
    sprintf(filename, "%d-%d-%d_%d.graph", R, S, NV, seed);

    for (igen = 0; igen <= ngen; igen++)
    {
        GA_advance(&ga, &ncross, &nmutation);

        if (igen % 1 == 0)
        {
            if (igen % 100 == 0) 
            {
                printf("#%8s %9s %9s %9s %9s %9s %9s %9s\n",
                       "gen", "emin", "favg", "fvar", "fmin", "fmax",
                       "ncross", "nmutation");
            }

            printf("%9d %9.0f %9.3g %9.3g %9.3g %9.3g %9d %9d\n",
                    igen,
                    ga.objective[ga.fittest],
                    ga.favg,
                    ga.fvar,
                    ga.fmin,
                    ga.fmax,
                    ncross,
                    nmutation
                  );

            fflush(stdout);
        }

        if (ga.objective[ga.fittest] < 1e-9) break;
    }

    decode(ga.chrom[ga.fittest], sp);
    R_save_graph(sp, filename);

    return EXIT_SUCCESS;
}
