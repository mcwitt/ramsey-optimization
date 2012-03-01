#include <stdio.h>
#include "ramsey.h"
#include "sga.h"

void decode(int chrom[], int sp[])
{
    int i;

    for (i = 0; i < NED; i++)
        sp[i] = (chrom[i] == 1) ? 1 : -1;
}

double SGA_objfunc(int chrom[])
{
    int sp[NED];

    decode(chrom, sp);
    return 1./(R_energy(sp)+1);
}

int main(int argc, char *argv[])
{
    SGA_indiv_t b1[SGA_MAXPOPSIZE];
    SGA_indiv_t b2[SGA_MAXPOPSIZE];
    SGA_indiv_t *newpop, *oldpop, *swap;
    SGA_params_t params;
    SGA_stats_t stats;
    char filename[256];
    int sp[NED];
    int igen, ngen, seed;

    if (argc != 6)
    {
        fprintf(stderr, "Usage: %s ngen npop pcross pmutate seed\n", argv[0]);
        fprintf(stderr, "Compiled for (%d, %d, %d)\n", R, S, NV);
        exit(EXIT_FAILURE);
    }

    ngen            = atoi(argv[1]);
    params.popsize  = atoi(argv[2]);
    params.pcross   = atof(argv[3]);
    params.pmutate  = atof(argv[4]);
    seed            = atoi(argv[5]);

    params.lchrom = NED;

    newpop = b1;
    oldpop = b2;

    stats.ncross = 0;
    stats.nmutation = 0;

    R_init(seed);
    SGA_init(seed);
    SGA_init_pop(oldpop, &params);

    sprintf(filename, "%d-%d-%d_%d.graph", R, S, NV, seed);
    printf("%9s %9s %9s %9s %9s %9s %9s\n",
            "gen", "avg", "var", "max", "min", "ncross", "nmutation");

    for (igen = 0; igen < ngen; igen++)
    {
        SGA_advance(oldpop, newpop, &params, &stats);

        printf("%9d %9.3g %9.3g %9.3g %9.3g %9d %9d\n",
                igen,
                stats.fitness_avg,
                stats.fitness_var,
                stats.fitness_max,
                stats.fitness_min,
                stats.ncross,
                stats.nmutation
              );
        fflush(stdout);

        if (stats.fitness_max > 0.999) break;

        /* swap pointers to population arrays so that oldpop becomes newpop
         * (former oldpop will be overwritten in next iteration) */
        swap = oldpop;
        oldpop = newpop;
        newpop = swap;
    }

    decode(stats.fittest->chrom, sp);
    R_save_graph(sp, filename);

    return EXIT_SUCCESS;
}
