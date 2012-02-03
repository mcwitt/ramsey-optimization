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
int emin;

int emax_demon_ini;
int nsweep_ini;
int nstage;
int nrun;
int itry, irun;
uint32_t seed;

void converge()
{
    int nsweep, nflip, nflip_sweep, emax_demon, e_demon, e_demon_av, istage,
        isweep, delta, j;

    nsweep = nsweep_ini;
    e_demon = emax_demon = emax_demon_ini;
    nflip_sweep = 0;

    for (istage = 0; istage < nstage; istage++)
    {
        nflip = 0;
        e_demon_av = 0;
        emax_demon = (int) ceil(emax_demon * ( 1. - istage/(nstage-1.) ));

        for (isweep = 0; isweep < nsweep; isweep++)
        {
            nflip_sweep = 0;

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
                    nflip_sweep++;
                }
            }

            e_demon_av += e_demon;
            if (! nflip_sweep) break;
            nflip += nflip_sweep;
        }

#ifdef FULL_OUTPUT
        /* print stats */
        printf("%3d %6d %6d %12d %12.2f %6.2f %6d\n",
                irun, itry, istage, emax_demon,
                (double) e_demon_av/(isweep+1),
                (double) nflip/NED/(isweep+1),
                (emin == INT_MAX) ? -1 : emin);
#endif
        if (! nflip_sweep) break;
    }
}

int main(int argc, char *argv[])
{
    char filename[256];
    int ntry[NRUN_MAX];
    int ntry_max;
    int imask;
    int converged;

    if (argc != 6 && argc != 7)
    {
        fprintf(stderr, "Usage: %s nsweep_ini emax_demon_ini nstage ntry_max"
                " seed [partial config]\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    nsweep_ini = atoi(argv[1]);
    emax_demon_ini = atoi(argv[2]);
    nstage = atoi(argv[3]);
    ntry_max = atoi(argv[4]);
    seed = atoi(argv[5]);

    sprintf(filename, "%d-%d-%d_%d.graph", R, S, NV, seed);
    
    R_init(seed);

    if (argc == 7)
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

    /* init "decades" */
    ntry[0] = 0; ntry[1] = 1; ntry[2] = 3;
    nrun = 3;
    while ((ntry[nrun] = 10*ntry[nrun-2]) <= ntry_max) nrun++;

    /* print header */
#ifdef FULL_OUTPUT
    printf("%3s %6s %6s %12s %12s %6s %6s\n",
            "run", "try", "stage", "emax_demon", "e_demon_av", "a.r.", "emin");
#else
    printf("%8s %8s\n", "N_try", "E_min");
#endif

    /* BEGIN SIMULATION */
    converged = 0;
    emin = INT_MAX;

    for (irun = 1; irun <= nrun; irun++)
    {
        for (itry = 0; itry < (ntry[irun] - ntry[irun-1]); itry++)
        {
            R_randomize(&r, imask);   /* randomize free spins */
            converge();

            if (r.en < emin)
            {
                emin = r.en;
                if (emin < WRITE_MAX) R_save_graph(r.sp, filename);
                if (emin == 0) { converged = 1; break; }
            }
        }

#ifndef FULL_OUTPUT
        printf("%8d %8d\n", ntry[irun], emin);
        fflush(stdout);
#endif
        if (converged) break;
    }

    R_finalize();
    return (emin == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
