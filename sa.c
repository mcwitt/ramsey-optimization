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
#include "dSFMT.h"
#include "ramsey.h"

#define WRITE_MAX 10  /* only save graph when energy is below this value */

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
    double T_max, T_min, T_mult, sweep_mult;
    int nsweep_min, nsweep_max, nstage, nrun;
    int irun, istage, isweep;
    int nsweep, nflip;
    int emin, emin_stage;
    int mask;
    uint32_t seed;

    if (argc != 8 && argc != 9)
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
    nrun = atoi(argv[6]);
    seed = atoi(argv[7]);

    T_mult = pow(T_min / T_max, 1./nstage);
    sweep_mult = pow((double) nsweep_max / nsweep_min, 1./nstage);
    sprintf(filename, "%d-%d-%d_%d.graph", R, S, NV, seed);
    
    dsfmt_init_gen_rand(&dsfmt, seed);
    R_init(seed);

    if (argc == 9)
    {
        /*
         * load configuration from file and set mask to prevent spins
         * specified in the input file from being randomized before each
         * iteration
         */

        mask = R_init_replica_from_file(&r, argv[8]);
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
        printf("# %3s %8s %8s %10s %8s %12s %8s\n",
            "run", "stage", "nsweep", "T", "a.r.", "emin_stage", "emin");

        R_randomize(&r, (double) R/(R+S), mask);   /* randomize free spins */
        nsweep = nsweep_min;
        T = T_max;

        for (istage = nstage; istage >= 0; istage--)
        {
            nflip = 0;
            emin_stage = INT_MAX;

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
            printf("%5d %8d %8d %10.2e %8.5f %12d %8d\n",
                    irun, istage, nsweep, T,
                    (double) nflip/NED/(isweep+1),
                    emin_stage, emin);
            fflush(stdout);

            if (emin == 0) break;
            T *= T_mult;
            nsweep *= sweep_mult;
        }

        if (emin == 0) break;
    }

    R_finalize();
    return (emin == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
