/*
 * Undefined constants
 *      NV  : number of vertices
 *      S   : clique size
 */

#define NT_MAX 16

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include "dSFMT.h"

#define M S-2
#define URAND() dsfmt_genrand_close_open(&rstate)
typedef unsigned long ULONG;

/* Replica-specific variables ************************************************/
typedef struct
{
    int s[NV][NV];  /* edge matrix (+1=blue, -1=red) */
    int h2[NV][NV]; /* field */
    int energy;
} rep_t;

/* Global variables **********************************************************/

int nsg;    /* number of subgraphs with S vertices (=binomial(N, S)) */
int nsg_fe; /* " involving a given edge (=binomial(N-2, S-2)) */

int *sub[NV][NV];   /* sub[i][j] is an array of subgraphs each containing S
                       vertices and including the edge (i, j) */

rep_t reps[NT_MAX]; /* storage for parallel tempering (PT) replicas */

/* pointers to PT replicas in order of increasing temperature */
rep_t *preps[NT_MAX];   

int nT;                 /* number of PT copies */
int nswaps[NT_MAX];     /* number of swaps between each pair of temperatures */
double T[NT_MAX];       /* array of temperatures */
double mbeta[NT_MAX];   /* negative inverse temperatures */

dsfmt_t rstate; /* state of random number generator (RNG) */
uint32_t rseed; /* seed used to initialize RNG */

/* binomial coefficient */
ULONG binomial(ULONG n, ULONG k)
{
    ULONG r = 1, d = n - k; 

    /* choose the smaller of k and n-k */
    if (d > k) { k = d; d = n - k; }
    while (n > k)
    {
        if (r >= ULONG_MAX / n) return 0;    /* overflown */
        r *= n--;
        /* divide as soon as possible to avoid overflow */
        while (d > 1 && !(r % d)) r /= d--;
    }
    return r;
}

void init_globals()
{
    int j, k;

    nsg = binomial(NV, S);
    nsg_fe = S*(S-1.)/(NV*(NV-1.))*nsg; /* (=binomial(NV-2, S-2)) */
    /* debug */
    /* assert(nsg_fe==binomial(NV-2,S-2)); */

    for (j = 0; j < NV; j++)
        for (k = 0; k < j; k++)
            sub[j][k] = (int*) malloc(nsg_fe * sizeof(int));

    /* fill in values */

}

void free_globals()
{
    int i, j;

    for (i = 0; i < NV; i++)
        for (j = 0; j < NV; j++)
            free(sub[i][j]);
}

/* initialize each replica with all spins +1 (i.e. all edges blue) */
void init_reps(rep_t *reps, rep_t **preps)
{
    rep_t *p;
    int iT, j, k;

    for (iT = 0; iT < nT; iT++)
    {
        p = &reps[iT];
        preps[iT] = p;
        p->energy = nsg;

        for (j = 0; j < NV; j++)
        {
            for (k = 0; k < j; k++)
            {
                p->s[j][k] = 1;
                p->h2[j][k] = nsg_fe;
            }
        }
    }   /* end of loop over temperatures */
}

/* flip a spin and update fields */
void flip(rep_t *p, int j, int k, int delta)
{
    int l, m, a[M];

    p->energy += delta;

    /* update local field at affected edges */
    if ((p->s[j][k] *= -1) == 1)
    {
    }
    else
    {
        /* same for combos of S vertices, subtract 1 */
    }
}

void sweep(rep_t **preps)
{
    rep_t *p;
    int iT, j, k, delta;

    for (iT = 0; iT < nT; iT++)
    {
        p = preps[iT];
        for (j = 0; j < NV; j++)
        {
            for (k = 0; k < j; k++)
            {
                /* compute energy difference of flip */
                delta = p->s[j][k]*p->h2[j][k];

                /* flip with Metropolis probability */
                if (delta < 0 || URAND() < exp(mbeta[iT]*delta))
                    flip(p, j, k, delta);
            }
        }
    }   /* end of loop over temperatures */
}

int main(int argc, char *argv[])
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
    {
        T[nT] = t;
        mbeta[nT] = 1./t;
        nT++;
    }

    init_globals();
    
    free_globals();

    return EXIT_SUCCESS;
}
