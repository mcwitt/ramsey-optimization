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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "dSFMT.h"
#include "qselect.h"
#include "ramsey.h"

#define WRITE_MAX 10  /* only save graph when energy is below this value */

#define RANDOM() dsfmt_genrand_close_open(&dsfmt)

R_replica_t r;
dsfmt_t dsfmt;

/* compute the discrete cumulative probabity distribution corresponding to
 * P(k) = 1/k^tau */
void set_cdf(double tau, double cdf[])
{
    double sum;
    int k;

    sum = 0.;
    for (k = 0; k < NED; k++)
        sum += (cdf[k] = pow(k+1, -tau));
    cdf[0] /= sum;
    for (k = 1; k < NED; k++)
        cdf[k] = cdf[k-1] + cdf[k]/sum;

    /*debug
    for (k = 0; k < NED; k++)
        printf("%f\n",cdf[k]);
    exit(0);*/
}

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
    double tau, cdf[NED];
    int lambda[NED];
    int nupdate, iupdate, j, k, emin;
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
    set_cdf(tau, cdf);

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

    printf("%12s %6s\n", "iter/N_edge", "emin");

    for (iupdate = 0; iupdate < nupdate; iupdate++)
    {

        /* compute fitness value for each spin */
        for (j = 0; j < NED; j++) lambda[j] = r.sp[j]*r.h2[j];

        /* choose random integer k in [0..NED-1] from distribution */
        k = bisect(cdf, RANDOM(), 0, NED);

        /* flip the kth "worst" spin */
        j = qselect_index(k, NED, lambda);
        r.en += lambda[j];
        R_flip(&r, j);

        if (r.en < emin)
        {
            emin = r.en;

            if (emin < WRITE_MAX)
            {
                R_save_graph(r.sp, filename);
                if (emin == 0) break;
            }
        }

        if (iupdate % NED == 0)
        {
            printf("%12d %6d\n", iupdate/NED, emin);
            fflush(stdout);
        }
    }

    printf("%5d %6d\n", iupdate, emin);
    R_finalize();
    return (emin == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
