/*
 * File: demon.c
 *
 * Author: Matt Wittmann <mwittman@ucsc.edu>
 *
 * Description: Annealed Demon Algorithm (ADA) Monte Carlo code which attempts
 * to minimize the number of r-cliques and s-independent sets of a graph given
 * that it has N_v vertices. Energy is defined to be the sum of the number of
 * r-cliques and s-independent sets. If for a given input (r, s, N_v) we find a
 * zero-energy state, this implies that R(r, s) > N_v.
 */

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ramsey.h"

#define WRITE_MAX 10  /* only save graph when energy is below this value */

R_replica_t r;
int e_demon;

int sweep(int emax_demon)
{
    int j, delta, nflip = 0;

    for (j = 0; j < NED; j++)
    {
        delta = r.sp[j]*r.h2[j];

        if (delta < e_demon)
        {
            r.en += delta;
            e_demon -= delta;
            R_flip(&r, j);
            if (e_demon > emax_demon) e_demon = emax_demon;
            nflip++;
        }
    }

    return nflip;
}

int main(int argc, char *argv[])
{
    char filename[256];
    double sweep_mult;
    int nsweep_min, nsweep_max, nstage, nrun;
    int irun, isweep;
    int nsweep, nflip, nflip_sweep = 0;
    int emax_demon, e_demon_av, emin, emin_stage;
    int mask;
    uint32_t seed;

    if (argc != 6 && argc != 7)
    {
        fprintf(stderr, "Usage: %s nsweep_min nsweep_max"
                " nstage nrun seed [partial_config]\n", argv[0]);
        fprintf(stderr, "Compiled for (%d, %d, %d)\n", R, S, NV);
        exit(EXIT_FAILURE);
    }

    nsweep_min = atoi(argv[1]);
    nsweep_max = atoi(argv[2]);
    nstage = atoi(argv[3]);
    nrun = atoi(argv[4]);
    seed = atoi(argv[5]);

    sweep_mult = pow((double) nsweep_max / nsweep_min, 1./nstage);
    sprintf(filename, "%d-%d-%d_%d.graph", R, S, NV, seed);
    
    R_init(seed);

    if (argc == 7)
    {
        /*
         * load configuration from file and set mask to prevent spins
         * specified in the input file from being randomized before each
         * iteration
         */

        mask = R_init_replica_from_file(&r, argv[6]);
        printf("Starting from configuration in %s. Mask = %d\n",
                filename, mask);
        assert(mask < NED);
    }
    else
    {
        R_init_replica(&r);
        mask = 0;
    }

    emin = INT_MAX;

    for (irun = 0; irun < nrun; irun++)
    {
        /* print column names */
        printf("# %3s %8s %12s %12s %8s %10s %12s %8s\n",
            "run", "nsweep", "emax_demon", "e_demon_av",
            "a.r.", "nflip/spin", "emin_stage", "emin");

        R_randomize(&r, (double) R/(R+S), mask);   /* randomize free spins */
        nsweep = nsweep_min;
        e_demon = emax_demon = nstage;

        while (emax_demon >= 0)
        {
            nflip = 0;
            e_demon_av = 0;
            emin_stage = INT_MAX;

            for (isweep = 0; isweep < nsweep; isweep++)
            {
                nflip_sweep = sweep(emax_demon);
                e_demon_av += e_demon;
                if (nflip_sweep == 0) break;
                nflip += nflip_sweep;

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
            printf("%5d %8d %12d %12.2f %8.5f %10.2f %12d %8d\n",
                    irun, nsweep, emax_demon,
                    (double) e_demon_av/(isweep+1),
                    (double) nflip/NED/(isweep+1),
                    (double) nflip/NED,
                    emin_stage, emin);
            fflush(stdout);

            if ((emin == 0) || (nflip_sweep == 0)) break;
            nsweep *= sweep_mult;
            emax_demon -= 1;
        }

        if (emin == 0) break;
    }

    R_finalize();
    return (emin == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
