/* 
 * File: ramsey.c
 *
 * Version 3
 *
 * Author: Matt Wittmann <mwittman@ucsc.edu>
 *
 * Part of ramsey3, a parallel tempering Monte Carlo code which attempts to
 * minimize the number of r-cliques and s-independent sets of a graph given
 * that it has N_v vertices. Energy is defined to be the sum of the number of
 * r-cliques and s-independent sets. If for a given input (r, s, N_v) we find a
 * zero-energy state, this implies that R(r, s) > N_v.
 *
 * Undefined constants (computed by compile.py)
 *      NV      : number of vertices
 *      S       : clique size
 *      NED     : number of edges (=NV(NV-1)/2)
 *      NSG     : number of subgraphs with S vertices (=binomial(NV, S))
 *      NSGFE   : number of subgraphs with S vertices including a given edge
 *                  (=binomial(NV-2, S-2))
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "ramsey.h"

/* Global variables **********************************************************/

int neds    = S*(S-1)/2;        /* number of edges in an S-subgraph */
int nedsm1  = S*(S-1)/2 - 1;
int nedsm2  = S*(S-1)/2 - 2;

/* sub[ei] lists the complete S-subgraphs that include edge ei */
/* edg[si] lists the edges of subgraph si */
int *sub[NED];
int *edg[NSG];

rep_t reps[MAX_NT]; /* storage for parallel tempering (PT) replicas */
int ri[MAX_NT];     /* replica indices in order of increasing temperature */

int nsweeps;    /* number of sweeps */
int min;        /* lowest energy found */

int nt;                 /* number of PT copies */
int nswaps[MAX_NT];     /* number of swaps between each pair of temperatures */
double T[MAX_NT];       /* array of temperatures */
double mbeta[MAX_NT];   /* negative inverse temperatures */

dsfmt_t rstate; /* state of random number generator (RNG) */
uint32_t rseed; /* seed used to initialize RNG */

#ifndef NOTIME
clock_t start;  /* start time */
#endif

void free_tabs()
{
    int i;

    for (i = 0; i < NED; i++)
        free(sub[i]);
    for (i = 0; i < NSG; i++)
        free(edg[i]);
}

int main(int argc, char *argv[])
{
    FILE *infile;
    double t;

    if (argc != 3 && argc != 4)
    {
        fprintf(stderr, "Usage: %s input_file seed [saved state]\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    /* init random number generator */
    rseed = atoi(argv[2]);
    dsfmt_init_gen_rand(&rstate, rseed);

    /* read temperatures from input file */
    nt = 0;
    infile = fopen(argv[1], "r");
    while (fscanf(infile, "%lf", &t) != EOF && nt < MAX_NT)
    {
        T[nt] = t;
        mbeta[nt] = -1./t;
        nt++;
    }
    assert(nt > 1);

    init_tabs();
    init_reps();

    if (argc == 4) load_state(argv[3]);

    run();
    
    free_tabs();

    return EXIT_SUCCESS;
}
