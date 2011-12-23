/*
 * Undefined constants
 *      NV: number of vertices
 *      NT_MAX: maximum number of temperatures
 */

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include "dSFMT.h"

/* Replica-specific variables ************************************************/
rep_t[NT_MAX] reps; /* storage for parallel tempering (PT) replicas */
typedef struct
{
    int am[N][N];   /* adjacency matrix */
    int energy;
} rep;

/* Global variables **********************************************************/
rep_t reps[NT_MAX]; /* storage for parallel tempering (PT) replicas */

/* pointers to PT replicas in order of increasing temperature */
rep_t *preps[NT_MAX];   

int nT;             /* number of PT copies */
double T[NT_MAX];   /* array of temperatures */
int nswaps[NT_MAX]; /* number of swaps between each adjacent pair */

dsfmt_t rstate;     /* state of random number generator */
uint32_t rseed;

void sweep(rep_t **preps)
{
    rep_t *p;
    double delta;
    int iT;

    for (iT = 0; iT < nT; iT++)
    {
        p = preps[iT];
    }
}

int main()
{
    FILE *infile;
    double t;

    if (argc != 3)
    {
        fprintf(stderr, "Usage: %s input_file seed\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    rseed = atoi(argv[2]);

    /* read temperatures from input file */
    nT = 0;
    infile = fopen(argv[1], "r");
    while (fscanf(infile, "%lf", &t) != EOF)
        T[nT++] = t;

    return EXIT_SUCCESS;
}
