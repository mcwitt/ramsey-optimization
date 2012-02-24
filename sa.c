/*
 * File: sa.c
 *
 * Author: Matt Wittmann <mwittman@ucsc.edu>
 *
 * Description: Simulated Annealing Monte Carlo code which attempts to minimize
 * the number of r-cliques and s-independent sets of a graph given that it has
 * N_v vertices. Energy is defined to be the sum of the number of r-cliques and
 * s-independent sets. If for a given input (r, s, N_v) we find a zero-energy
 * state, this implies that R(r, s) > N_v.
 */

#include <limits.h>
#include <math.h>
#include "fastexp.h"
#include "ramsey.h"

#define WRITE_MAX 10  /* only save graph when energy is below this value */

rep_t r;
double T;

int sweep()
{
    int j, delta, nflip = 0;

    for (j = 0; j < NED; j++)
    {
        delta = r.sp[j]*r.h2[j];

        if (delta <= 0 || R_RAND() < EXP(-delta/T))
        {
            R_flip(&r, j);
            r.en += delta;
            nflip++;
        }
    }

    return nflip;
}

int main(int argc, char *argv[])
{
    char filename[256];
    double T_ini, sweep_mult;
    int nsweep_ini, nstage, nrun;
    int irun, istage, isweep;
    int nsweep, nflip;
    int emin, emin_stage;
    int mask;
    uint32_t seed;

    if (argc != 7 && argc != 8)
    {
        fprintf(stderr, "Usage: %s T_ini nsweep_ini sweep_mult"
                " nstage nrun seed [partial_config]\n", argv[0]);
        fprintf(stderr, "Compiled for (%d, %d, %d)\n", R, S, NV);
        exit(EXIT_FAILURE);
    }

    T_ini = atof(argv[1]);
    nsweep_ini = atoi(argv[2]);
    sweep_mult = atof(argv[3]);
    nstage = atoi(argv[4]);
    nrun = atoi(argv[5]);
    seed = atoi(argv[6]);

    sprintf(filename, "%d-%d-%d_%d.graph", R, S, NV, seed);
    
    R_init(seed);

    if (argc == 8)
    {
        /*
         * load configuration from file and set mask to prevent spins
         * specified in the input file from being randomized before each
         * iteration
         */

        mask = R_init_replica_from_file(&r, argv[7]);
        printf("Starting from configuration in %s. Mask = %d\n",
                filename, mask);
        assert(mask < NED);
    }
    else
    {
        R_init_replica(&r);
        mask = 0;
    }

    /* BEGIN SIMULATION */
    emin = INT_MAX;

    for (irun = 0; irun < nrun; irun++)
    {
        /* print column names */
        printf("# %3s %8s %8s %8s %8s %12s %8s\n",
            "run", "stage", "nsweep", "T", "a.r.", "emin_stage", "emin");

        R_randomize(&r, (double) R/(R+S), mask);   /* randomize free spins */
        nsweep = nsweep_ini;
        T = T_ini;

        for (istage = nstage; istage >= 0; istage--)
        {
            nflip = 0;
            emin_stage = INT_MAX;
            T *= (double) istage / nstage;

            for (isweep = 0; isweep < nsweep; isweep++)
            {
                nflip += sweep();

                if (r.en < emin_stage)
                {
                    emin_stage = r.en;

                    if (emin_stage < emin)
                    {
                        emin = emin_stage;

                        if (emin < WRITE_MAX)
                        {
                            R_save_graph(r.sp, filename);
                            if (emin == 0) break;
                        }
                    }
                }
            }

            /* print stage stats */
            printf("%5d %8d %8d %8.5f %8.5f %12d %8d\n",
                    irun, istage, nsweep, T,
                    (double) nflip/NED/(isweep+1),
                    emin_stage, emin);
            fflush(stdout);

            if (emin == 0) break;
            nsweep *= sweep_mult;
        }

        if (emin == 0) break;
    }

    R_finalize();
    return (emin == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
