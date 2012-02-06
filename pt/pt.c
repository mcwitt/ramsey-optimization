/*
 * File: pt.c
 *
 * Author: Matt Wittmann <mwittman@ucsc.edu>
 *
 * Description: Parallel tempering Monte Carlo code which attempts to minimize
 * the number of r-cliques and s-independent sets of a graph given that it has
 * N_v vertices. Energy is defined to be the sum of the number of r-cliques and
 * s-independent sets. If for a given input (r, s, N_v) we find a zero-energy
 * state, this implies that R(r, s) > N_v.
 *
 * To compile, run python compile.py.
 */

#define MAX_NT      32  /* maximum number of parallel tempering replicas */
#define WRITE_MAX   10  /* only save graph when energy is below this value */

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "ramsey.h"

rep_t reps[MAX_NT]; /* storage for parallel tempering (PT) replicas */
int ri[MAX_NT];     /* replica indices in order of increasing temperature */

int nsweeps;        /* number of sweeps */
int emin;           /* lowest energy found */
int max_sweeps;     /* number of sweeps to do before giving up */

int nT;                 /* number of PT copies */
int nswaps[MAX_NT];     /* number of swaps between each pair of temperatures */
double T[MAX_NT];       /* array of temperatures */
double mbeta[MAX_NT];   /* negative inverse temperatures */

uint32_t seed; /* seed used to initialize RNG */

#ifndef NOTIME
clock_t start;  /* start time */
#endif

/* Update each spin once */
void sweep()
{
    rep_t *p;
    int iT, j, delta;

    for (iT = 0; iT < nT; iT++)
    {
        p = &reps[ri[iT]];

        for (j = 0; j < NED; j++)
        {
            /* compute energy difference of flip */
            delta = p->h2[j]*p->sp[j];

            /* flip with Metropolis probability */
            if (delta <= 0 || R_RAND() < exp(mbeta[iT]*delta))
            {
                p->sp[j] *= -1;
                p->en += delta;
                R_update_fields(p, j);
            }
        }
    }
}

/* Attempt parallel tempering swaps */
void temper()
{
    double logar;
    int iT, copy;

    for (iT = 1; iT < nT; iT++)
    {
        logar = (reps[ri[iT-1]].en - reps[ri[iT]].en)
            * (mbeta[iT] - mbeta[iT-1]);

        if (R_RAND() < exp(logar))
        {
            /* do PT swap */
            copy = ri[iT-1];
            ri[iT-1] = ri[iT];
            ri[iT] = copy;
            nswaps[iT]++;
        }
    }
}

void print_status()
{
    int iT;
#ifndef NOTIME
    int trun;
#endif

    printf("\n");
    printf("min. energy     : %d\n", emin);
    printf("# of sweeps     : %d\n", nsweeps);
#ifndef NOTIME
    trun = (clock() - start)/CLOCKS_PER_SEC;
    printf("time running    : %d seconds\n", trun);
    printf("sweep rate      : %.2f / s\n", (float) nsweeps/trun);
#endif
    for (iT = 0; iT < nT; iT++)
        printf("%5d ", reps[ri[iT]].en);
    printf("\n");
    for (iT = 0; iT < nT; iT++)
        printf("%5.2f ", T[iT]);
    printf("\n   ");
    for (iT = 1; iT < nT; iT++)
        printf("%5.2f ", (float) nswaps[iT]/nsweeps);
    printf("\n");
    fflush(stdout);
}

int main(int argc, char *argv[])
{
    FILE *infile;
    rep_t *p;
    char filename[256];
    double t;
    int iT;

#ifndef NOTIME
    start = clock();
#endif

    if (argc != 4 && argc != 5)
    {
        fprintf(stderr, "Usage: %s T_file max_sweeps"
               " seed [initial config]\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    max_sweeps = atoi(argv[2]);
    seed = atoi(argv[3]);

    /* READ TEMPERATURES */
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


    /* INITIALIZE SIMULATION */
    R_init(seed);

    /* INITIALIZE REPLICAS */
    if (argc == 5) /* initial configuration specified */
    {
        R_init_replica_from_file(reps, argv[4]);
        for (iT = 0; iT < nT; iT++)
            reps[iT] = reps[0];
    }
    else
        for (iT = 0; iT < nT; iT++)
            R_init_replica_random(&reps[iT]);

    for (iT = 0; iT < nT; iT++)
    {
        ri[iT] = iT;
        nswaps[iT] = 0;
    }

    /* BEGIN SIMULATION */
    sprintf(filename, "%d-%d-%d_%d.graph", R, S, NV, seed);
    emin = INT_MAX;
    nsweeps = 0;

    while ((nsweeps < max_sweeps) && (emin > 0))
    {
        sweep();
        nsweeps++;

        temper();

        for (iT = 0; iT < nT; iT++)
        {
            p = &reps[ri[iT]];
            if (p->en < emin)
            {
                emin = p->en;
                print_status();

                if (emin < WRITE_MAX)
                {
                    R_save_graph(p->sp, filename);
                    if (emin == 0) break;
                }
            }
        }
    }

    print_status();
    R_finalize();

    return (emin == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
