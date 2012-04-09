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
    char filename[64];
    double sweep_mult, nflip = 0.;
    int nsweep_min, nsweep_max, nsweep, nstage, isweep, emax_demon, mask,
        nflip_sweep = 0;
    uint32_t seed;

    if (argc != 5 && argc != 6)
    {
        fprintf(stderr, "Usage: %s nsweep_min nsweep_max"
                " nstage seed [partial_config]\n", argv[0]);
        fprintf(stderr, "Compiled for (%d, %d, %d)\n", R, S, NV);
        exit(EXIT_FAILURE);
    }

    nsweep_min = atoi(argv[1]);
    nsweep_max = atoi(argv[2]);
    nstage = atoi(argv[3]);
    seed = atoi(argv[4]);

    sweep_mult = pow((double) nsweep_max / nsweep_min, 1./nstage);
    sprintf(filename, "%d-%d-%d_%d.graph", R, S, NV, seed);
    
    R_init(seed);

    if (argc == 6)
    {
        /*
         * load configuration from file and set mask to prevent spins
         * specified in the input file from being randomized before each
         * iteration
         */

        mask = R_init_replica_from_file(&r, argv[5]);
        fprintf(stderr, "Starting from configuration in %s. Mask = %d\n",
                filename, mask);
        assert(mask < NED);
    }
    else
    {
        R_init_replica(&r);
        mask = 0;
    }

    while (r.en > 0)
    {
        R_randomize(&r, (double) R/(R+S), mask);   /* randomize free spins */
        nsweep = nsweep_min;
        e_demon = emax_demon = nstage;

        while (emax_demon >= 0)
        {
            for (isweep = 0; isweep < nsweep; isweep++)
            {
                nflip_sweep = sweep(emax_demon);
                if (nflip_sweep == 0) break;
                nflip += (double) nflip_sweep;
                if (r.en == 0) break;
            }

            if ((r.en == 0) || (nflip_sweep == 0)) break;
            nsweep *= sweep_mult;
            emax_demon -= 1;
        }
    }

    printf("%g\n", nflip/NED);
    R_save_graph(r.sp, filename);
    R_finalize();
    return EXIT_SUCCESS;
}
