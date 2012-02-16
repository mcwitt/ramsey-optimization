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
#include "ramsey2.h"

#define WRITE_MAX 500.  /* only save graph when energy is below this value */
#define NRUN_MAX 100

rep_t r;
double e_demon;

#if defined(LINEAR)
double er_ini[] = {1.00, 0.83, 0.67, 0.50, 0.33, 0.17, 0.00};
double es_ini[] = {
    1.00, 0.93, 0.87, 0.80, 0.73, 0.67, 0.60, 0.53, 0.47, 0.40, 0.33,
    0.27, 0.20, 0.13, 0.07, 0.00
};
#elif defined(QUADRATIC)
double er_ini[] = {1.00, 0.69, 0.44, 0.25, 0.11, 0.03, 0.00};
double es_ini[] = {
    1.00, 0.87, 0.75, 0.64, 0.54, 0.44, 0.36, 0.28, 0.22, 0.16, 0.11, 0.07,
    0.04, 0.02, 0.00, 0.00 
};
#else
double er_ini[] = {1., 0., 0., 0., 0., 0., 0.};
double es_ini[] = {
    1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
    0., 0., 0., 0., 0., 0.
};
#endif

void print_header()
{
    printf("# %3s %8s %8s %12s %12s %8s %10s %12s %8s\n",
            "run", "stage", "nsweep", "emax_demon", "e_demon_av",
            "a.r.", "nflip/spin", "emin_stage", "emin");
}

int sweep(int emax_demon)
{
    int j, nflip = 0;
    double delta;

    for (j = 0; j < NED; j++)
    {
        delta = R_flip_energy(&r, j);

        if (delta < e_demon)
        {
            r.en += delta;
            e_demon -= delta;
            r.sp[j] *= -1;
            R_update(&r, j);
            if (e_demon > emax_demon) e_demon = emax_demon;
            nflip += 1;
        }
    }

    return nflip;
}

int main(int argc, char *argv[])
{
    char filename[256];
    double er[NEDR+1], es[NEDS+1];
    int nsweep_ini, nstage, nrun;
    double emax_demon_ini, sweep_mult;

    int imask, irun, istage, isweep, j;
    int nsweep, nflip, nflip_sweep = 0;
    double emax_demon, e_demon_av, emin, emin_stage, e_mult;
    uint32_t seed;

    if (argc != 7 && argc != 8)
    {
        fprintf(stderr, "Usage: %s emax_demon_ini nsweep_ini sweep_mult"
                " nstage nrun seed [partial_config]\n", argv[0]);
        fprintf(stderr, "Compiled for (%d, %d, %d)\n", R, S, NV);
        exit(EXIT_FAILURE);
    }

    emax_demon_ini = atof(argv[1]);
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

        imask = R_init_replica_from_file(&r, argv[7]);
        assert(imask > 0);
    }
    else
    {
        R_init_replica(&r);
        imask = NED;
    }

    /* BEGIN SIMULATION */
    emin = 10e9;

    for (irun = 0; irun < nrun; irun++)
    {
        print_header();
        R_randomize(&r, (double) R/(R+S), imask);   /* randomize free spins */
        nsweep = nsweep_ini;
        e_demon = emax_demon = emax_demon_ini;

        for (j = 0; j < NEDR+1; j++) er[j] = er_ini[j];
        for (j = 0; j < NEDS+1; j++) es[j] = es_ini[j];
        R_set_energies(&r, er, es);

        for (istage = 0; istage < nstage; istage++)
        {
            nflip = 0;
            e_demon_av = 0.;
            emin_stage = 10e9;
            e_mult = 1. - istage/(nstage-1.);
#ifdef DEBUG
            emax_demon = (int) (emax_demon * e_mult);
#else
            emax_demon *= e_mult;
            for (j = 1; j < NEDR+1; j++) er[j] *= e_mult;
            for (j = 1; j < NEDS+1; j++) es[j] *= e_mult;
            R_set_energies(&r, er, es);
#endif

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
                            if (emin < 1.) break;
                        }
                    }
                }
            }

            /* print stats */
            printf("%5d %8d %8d %12.2f %12.2f %8.5f %10.2f %12.2f %8.2f\n",
                    irun, istage, isweep, emax_demon,
                    (double) e_demon_av/(isweep+1),
                    (double) nflip/NED/(isweep+1),
                    (double) nflip/NED,
                    emin_stage, emin);
            fflush(stdout);

            if ((emin < 1.) || (nflip_sweep == 0)) break;

            nsweep *= sweep_mult;
        }

        if (emin < 1.) break;
    }

    R_finalize();
    return (emin < 1.) ? EXIT_SUCCESS : EXIT_FAILURE;
}
