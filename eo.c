/*
 * File: eo.c
 *
 * Author: Matt Wittmann <mwittman@ucsc.edu>
 *
 * Description: Extremal Optimization (EO) Monte Carlo code which attempts
 * to minimize the number of r-cliques and s-independent sets of a graph given
 * that it has N_v vertices. Energy is defined to be the sum of the number of
 * r-cliques and s-independent sets. If for a given input (r, s, N_v) we find a
 * zero-energy state, this implies that R(r, s) > N_v.
 */

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "ramsey.h"

#define WRITE_MAX 10  /* only save graph when energy is below this value */

#define RANDOM() dsfmt_genrand_close_open(&dsfmt)

R_replica_t r;

/* return the insertion point for x to maintain sorted order of a */
int bisect(double a[], double x, int l, int r)
{
    int mid;

    while (r-l > 1)
    {
        mid = (l+r)/2;
        if (a[mid] < x) l = mid;
        else r = mid;
    }

    return (x < a[l]) ? l : r;
}

int main(int argc, char *argv[])
{
    char filename[256];
    double tau;
    int nupdate, iupdate, emin;
    uint32_t seed;

    if (argc != 4 && argc != 5)
    {
        fprintf(stderr, "Usage: %s tau nupdate seed "
                " [partial_config]\n", argv[0]);
        fprintf(stderr, "Compiled for (%d, %d, %d)\n", R, S, NV);
        exit(EXIT_FAILURE);
    }

    tau = atof(argv[1]);
    nupdate = atoi(argv[2]);
    seed = atoi(argv[3]);

    sprintf(filename, "%d-%d-%d_%d.graph", R, S, NV, seed);
    
    R_init(seed);
    dsfmt_init_gen_rand(&dsfmt, seed);

    if (argc == 5)
    {   /* load configuration from file */
        R_init_replica_from_file(&r, argv[4]);
        printf("Starting from configuration in %s\n", filename);
    }
    else
    {
        R_init_replica(&r);
        R_randomize(&r, (double) R/(R+S), 0);
    }

    emin = INT_MAX;

    /* 
    printf("# %3s %8s %12s %12s %8s %10s %12s %8s\n",
        "run", "nsweep", "emax_demon", "e_demon_av",
        "a.r.", "nflip/spin", "emin_stage", "emin");
    */

    for (iupdate = 0; iupdate < nupdate; iupdate++)
    {
        

        /*
        printf("%5d %8d %12d %12.2f %8.5f %10.2f %12d %8d\n",
                irun, nsweep, emax_demon,
                (double) e_demon_av/(isweep+1),
                (double) nflip/NED/(isweep+1),
                (double) nflip/NED,
                emin_stage, emin);
        fflush(stdout);
        */


        if (emin == 0) break;
    }

    R_finalize();
    return (emin == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
