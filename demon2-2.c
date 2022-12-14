/*
 * File: demon2.c
 *
 * Author: Matt Wittmann <mwittman@ucsc.edu>
 *
 * Description: Soft-energy version of demon.c. In this version the potential
 * stage is advanced only when the system becomes trapped in a local minimum.
 */

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ramsey2.h"

#define WRITE_MAX 500.  /* only save graph when energy is below this value */
#define NRUN_MAX 100

R_replica_t r;
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
#warning Potential shape not defined. Defaulting to square.
double er_ini[] = {1., 0., 0., 0., 0., 0., 0.};
double es_ini[] = {
    1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
    0., 0., 0., 0., 0., 0.
};
#endif

void print_header()
{
    printf("# %3s %8s %8s %8s %8s %12s %12s %8s %10s %12s %12s\n",
            "run", "try", "vstage", "dstage", "nsweep", "emax_demon", "e_demon_av",
            "a.r.", "nflip/spin", "emin_try", "emin");
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
    char filename[64];
    double er[NEDR+1], es[NEDS+1];
    int nsweep_ini, nrun, nvstage, ndstage;
    double emax_demon_ini, sweep_mult;

    int mask;
    int irun, itry, isweep, j;
    int vstage, dstage;
    int nsweep, nflip, nflip_sweep = 0;
    double emax_demon, e_demon_av, emin, emin_try;
    uint32_t seed;

    if (argc != 8 && argc != 9)
    {
        fprintf(stderr, "Usage: %s emax_demon_ini nsweep_ini sweep_mult"
                " nvstage ndstage nrun seed [partial_config]\n", argv[0]);
        fprintf(stderr, "Compiled for (%d, %d, %d)\n", R, S, NV);
        exit(EXIT_FAILURE);
    }

    emax_demon_ini = atof(argv[1]);
    nsweep_ini = atoi(argv[2]);
    sweep_mult = atof(argv[3]);
    nvstage = atoi(argv[4]);
    ndstage = atoi(argv[5]);
    nrun = atoi(argv[6]);
    seed = atoi(argv[7]);

    sprintf(filename, "%d-%d-%d_%d.graph", R, S, NV, seed);
    
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

    /* BEGIN SIMULATION */
    emin = 10e9;

    for (irun = 0; irun < nrun; irun++)
    {
        print_header();
        R_randomize(&r, (double) R/(R+S), mask);   /* randomize free spins */
        nsweep = nsweep_ini;
        e_demon = emax_demon = emax_demon_ini;
        vstage = nvstage;
        dstage = ndstage;
        for (j = 0; j < NEDR+1; j++) er[j] = er_ini[j];
        for (j = 0; j < NEDS+1; j++) es[j] = es_ini[j];
        R_set_energies(&r, er, es);

        for (itry = 0; ; itry++)
        {
            nflip = 0;
            e_demon_av = 0.;
            emin_try = 10e9;

            for (isweep = 0; isweep < nsweep; isweep++)
            {
                nflip_sweep = sweep(emax_demon);
                e_demon_av += e_demon;

                if (r.en < emin_try)
                {
                    emin_try = r.en;

                    if (emin_try < emin)
                    {
                        emin = emin_try;

                        if (emin < WRITE_MAX)
                        {
                            R_save_graph(r.sp, filename);
                            if (emin < 1.) break;
                        }
                    }
                }

                if (nflip_sweep == 0) break;
                nflip += nflip_sweep;
            }

            /* print stats */
            printf("%5d %8d %8d %8d %8d %12.2f %12.2f %8.5f %10.2f %12.2f %12.2f\n",
                    irun, itry, vstage, dstage, isweep, emax_demon,
                    (double) e_demon_av/(isweep+1),
                    (double) nflip/NED/(isweep+1),
                    (double) nflip/NED,
                    emin_try, emin);
            fflush(stdout);

            if (emin < 1.) break;

            if (nflip_sweep == 0)
            {
                if ((--vstage) < 0) break;
                for (j = 1; j < NEDR+1; j++) er[j] *= (double) vstage/nvstage;
                for (j = 1; j < NEDS+1; j++) es[j] *= (double) vstage/nvstage;
                R_set_energies(&r, er, es);
            }
            else
            {
                if ((--dstage) < 0) break;
                emax_demon *= (double) dstage/ndstage;
                nsweep *= sweep_mult;
            }
        }

        if (emin < 1.) break;
    }

    R_finalize();
    return (emin < 1.) ? EXIT_SUCCESS : EXIT_FAILURE;
}
