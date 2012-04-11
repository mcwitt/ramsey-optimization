/*
 * File: pt_s.c
 *
 * Author: Matt Wittmann <mwittman@ucsc.edu>
 *
 * Description: Parallel Tempering (PT) Monte Carlo code which attempts to
 * minimize the number of r-cliques and s-independent sets of a graph given
 * that it has N_v vertices. Energy is defined to be the sum of the number of
 * r-cliques and s-independent sets. If for a given input (r, s, N_v) we find a
 * zero-energy state, this implies that R(r, s) > N_v.
 */

#define MAX_NT      32  /* maximum number of parallel tempering replicas */

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "dSFMT.h"
#include "ramsey.h"

#define RANDOM() dsfmt_genrand_close_open(&dsfmt)

R_replica_t reps[MAX_NT];   /* storage for parallel tempering (PT) replicas */
int ri[MAX_NT];             /* replica indices in order of increasing temperature */

int emin;               /* lowest energy found */
int nT;                 /* number of PT copies */
double T[MAX_NT];       /* array of temperatures */
double mbeta[MAX_NT];   /* negative inverse temperatures */

dsfmt_t dsfmt;
uint32_t seed; /* seed used to initialize RNG */

/* sweep all replicas, return number of flips */
int sweep()
{
    R_replica_t *p;
    int iT, j, delta, nflip = 0;

    for (iT = 0; iT < nT; iT++)
    {
        p = &reps[ri[iT]];

        for (j = 0; j < NED; j++)
        {
            /* compute energy difference of flip */
            delta = p->h2[j]*p->sp[j];

            /* flip with Metropolis probability */
            if (delta <= 0 || RANDOM() < exp(mbeta[iT]*delta))
            {
                p->en += delta;
                R_flip(p, j);
                nflip++;
            }
        }
    }

    return nflip;
}

/* attempt parallel tempering swaps */
void temper()
{
    double logar;
    int iT, copy;

    for (iT = 1; iT < nT; iT++)
    {
        logar = (reps[ri[iT-1]].en - reps[ri[iT]].en)
            * (mbeta[iT] - mbeta[iT-1]);

        if (RANDOM() < exp(logar))
        {
            /* do PT swap */
            copy = ri[iT-1];
            ri[iT-1] = ri[iT];
            ri[iT] = copy;
        }
    }
}

int main(int argc, char *argv[])
{
    FILE *infile;
    R_replica_t *p;
    char filename[64];
    double t, nflip = 0.;
    int iT;

    if (argc != 3 && argc != 4)
    {
        fprintf(stderr, "Usage: %s T_file seed [initial config]\n", argv[0]);
        fprintf(stderr, "Compiled for (%d, %d, %d)\n", R, S, NV);
        exit(EXIT_FAILURE);
    }

    seed = atoi(argv[2]);

    /* read temperatures */
    if (! (infile = fopen(argv[1], "r")))
    {
        fprintf(stderr, "Error opening file: %s\n", argv[1]);
        exit(EXIT_FAILURE);
    }

    nT = 0;
    while (fscanf(infile, "%lf", &t) != EOF && nT < MAX_NT)
    {
        T[nT] = t;
        mbeta[nT] = -1./t;
        nT++;
    }
    assert(nT > 1);


    /* initialize simulation */
    dsfmt_init_gen_rand(&dsfmt, seed);
    R_init(seed);

    /* initialize replicas */
    if (argc == 4) /* initial configuration specified */
    {
        R_init_replica_from_file(reps, argv[3]);
        for (iT = 0; iT < nT; iT++) reps[iT] = reps[0];
    }
    else for (iT = 0; iT < nT; iT++)
    {
        R_init_replica(&reps[iT]);
        R_randomize(&reps[iT], (double) R/(R+S), 0);
    }

    for (iT = 0; iT < nT; iT++) ri[iT] = iT;

    /* begin simulation */
    sprintf(filename, "%d-%d-%d_%d.graph", R, S, NV, seed);
    emin = INT_MAX;

    while (emin > 0)
    {
        nflip += (double) sweep();
        temper();

        for (iT = 0; iT < nT; iT++)
        {
            p = &reps[ri[iT]];

            if (p->en < emin)
            {
                emin = p->en;

                if (emin == 0)
                {
                    R_save_graph(p->sp, filename);
                    break;
                }
            }
        }
    }

    printf("%g\n", nflip/NED);
    R_finalize();

    return EXIT_SUCCESS;
}
