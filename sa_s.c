/*
 * File: sa_s.c
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
#include "dSFMT.h"
#include "ramsey.h"

/*#include "fastexp.h"*/
#define EXP(x) exp(x)
#define RANDOM() dsfmt_genrand_close_open(&dsfmt)

R_replica_t r;
double T;

dsfmt_t dsfmt;

int sweep()
{
    int j, delta, nflip = 0;

    for (j = 0; j < NED; j++)
    {
        delta = r.sp[j]*r.h2[j];

        if (delta <= 0 || RANDOM() < EXP(-delta/T))
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
    char filename[64];
    double T_max, T_min, T_mult, sweep_mult, nflip = 0.;
    int nsweep_min, nsweep_max, nsweep, nstage, istage, isweep, mask;
    uint32_t seed;

    if (argc != 7 && argc != 8)
    {
        fprintf(stderr, "Usage: %s T_max T_min nsweep_min nsweep_max"
                " nstage nrun seed [partial_config]\n", argv[0]);
        fprintf(stderr, "Compiled for (%d, %d, %d)\n", R, S, NV);
        exit(EXIT_FAILURE);
    }

    T_max = atof(argv[1]);
    T_min = atof(argv[2]);
    nsweep_min = atoi(argv[3]);
    nsweep_max = atoi(argv[4]);
    nstage = atoi(argv[5]);
    seed = atoi(argv[6]);

    T_mult = pow(T_min / T_max, 1./nstage);
    sweep_mult = pow((double) nsweep_max / nsweep_min, 1./nstage);
    sprintf(filename, "%d-%d-%d_%d.graph", R, S, NV, seed);
    
    dsfmt_init_gen_rand(&dsfmt, seed);
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

    while (r.en > 0)
    {
        R_randomize(&r, (double) R/(R+S), mask);   /* randomize free spins */
        nsweep = nsweep_min;
        T = T_max;

        for (istage = nstage; istage >= 0; istage--)
        {
            for (isweep = 0; isweep < nsweep; isweep++)
            {
                nflip += (double) sweep();
                if (r.en == 0) break;
            }

            if (r.en == 0) break;
            T *= T_mult;
            nsweep *= sweep_mult;
        }
    }

    printf("%g\n", nflip/NED);
    R_save_graph(r.sp, filename);
    R_finalize();

    return EXIT_SUCCESS;
}
