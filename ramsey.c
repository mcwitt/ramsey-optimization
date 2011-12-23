/*
 * Undefined constants
 *      NV  : number of vertices
 *      R   : blue clique size
 *      S   : red clique size
 */

#define NT_MAX 16

#define URAND() dsfmt_genrand_close_open(&dsfmt)

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <limits.h>
#include "dSFMT.h"

typedef unsigned long ULONG;

/* Replica-specific variables ************************************************/
typedef struct
{
    int s[NV][NV];   /* adjacency matrix */
    int energy;
} rep_t;

/* Global variables **********************************************************/
rep_t reps[NT_MAX]; /* storage for parallel tempering (PT) replicas */

/* pointers to PT replicas in order of increasing temperature */
rep_t *preps[NT_MAX];   

int nT;             /* number of PT copies */
double T[NT_MAX];   /* array of temperatures */
int nswaps[NT_MAX]; /* number of swaps between each adjacent pair */

dsfmt_t rstate;     /* state of random number generator */
uint32_t rseed;

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

void init_reps(rep_t *reps, rep_t **preps)
{
    rep_t *p;
    int iT, j, k, e0, h0;

    /* initialize each replica in a simple state with all edges blue */

    /* local "field" at each edge in initial state */
    h20 = binomial(NV-2, R-2);

    /* total energy in initial state */
    e0 = NV*(NV-1)/(R*(R-1))*h20;

    for (iT = 0; iT < nT; iT++)
    {
        p = &reps[iT];
        preps[iT] = p;
        p->energy = e0;

        for (j = 0; j < NV; j++) {
            for (k = 0; k < j; k++) {
                p->s[j][k] = 1;
                p->h2[j][k] = h20;
            }
        }
    }   /* end of loop over temperatures */
}

/* flip a spin and update fields */
void flip(rep_t *p, int j, int k, int delta)
{
    p->energy += delta;

    /* update local field at affected edges */
    if ((p->s[j][k] *= -1) == 1)
    {
        /* DO UPDATE */
        /* loop over combos of R vertices */
        /* if a clique except for jk, add 1 */
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
                delta = p->s[j][k]*h2[j][k];

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
        T[nT++] = t;

    return EXIT_SUCCESS;
}
