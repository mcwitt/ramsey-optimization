/*
 * File: demon.c
 *
 * Author: Matt Wittmann <mwittman@ucsc.edu>
 *
 * Description:
 *
 * To compile, run python compile.py.
 */

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include "ramsey.h"

#define WRITE_MAX 10  /* only save graph when energy is below this value */
#define NRUN_MAX 100

rep_t r;
int e_demon;

void print_header()
{
    printf("# %3s %8s %8s %12s %12s %8s %10s %12s %8s\n",
            "run", "stage", "nsweep", "emax_demon", "e_demon_av",
            "a.r.", "nflip/spin", "emin_stage", "emin");
}

void sweep(int emax_demon, int *nflip)
{
    int j, delta;

    *nflip = 0;

    for (j = 0; j < NED; j++)
    {
        delta = r.sp[j]*r.h2[j];

        if (delta < e_demon)
        {
            r.en += delta;
            e_demon -= delta;
            r.sp[j] *= -1;
            R_update_fields(&r, j);
            if (e_demon > emax_demon) e_demon = emax_demon;
            *nflip += 1;
        }
    }
}

int main(int argc, char *argv[])
{
    char filename[256];
    int emax_demon_ini, nsweep_ini, nstage, nrun;
    int imask, irun, istage, isweep;
    int nsweep, nflip, nflip_sweep = 0;
    int emax_demon, e_demon_av, emin, emin_stage;
    double sweep_mult;
    uint32_t seed;

    if (argc != 7 && argc != 8)
    {
        fprintf(stderr, "Usage: %s emax_demon_ini nsweep_ini sweep_mult"
                " nstage nrun seed [partial config]\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    emax_demon_ini = atoi(argv[1]);
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
         * load configuration from file and set imask to prevent spins
         * specified in the input file from being randomized before each
         * iteration
         */

        imask = R_init_replica_from_file(&r, argv[4]);
        assert(imask > 0);
    }
    else
    {
        R_init_replica(&r);
        imask = NED;
    }

    /* BEGIN SIMULATION */
    emin = INT_MAX;

    for (irun = 0; irun < nrun; irun++)
    {
        print_header();
        R_randomize(&r, imask);   /* randomize free spins */
        nsweep = nsweep_ini;
        e_demon = emax_demon = emax_demon_ini;

        for (istage = 0; istage < nstage; istage++)
        {
            nflip = 0;
            e_demon_av = 0;
            emin_stage = INT_MAX;
            emax_demon = (int) (emax_demon * ( 1. - istage/(nstage-1.) ));

            for (isweep = 0; isweep < nsweep; isweep++)
            {
                sweep(emax_demon, &nflip_sweep);
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

            /* print stats */
            printf("%5d %8d %8d %12d %12.2f %8.5f %10.2f %12d %8d\n",
                    irun, istage, nsweep, emax_demon,
                    (double) e_demon_av/(isweep+1),
                    (double) nflip/NED/(isweep+1),
                    (double) nflip/NED,
                    emin_stage, emin);
            fflush(stdout);

            if ((emin == 0) || (nflip_sweep == 0)) break;

            nsweep *= sweep_mult;
        }

        if (emin == 0) break;
    }

    R_finalize();
    return (emin == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
