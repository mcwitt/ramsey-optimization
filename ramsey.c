/*
 * Undefined constants
 *      NV  : number of vertices
 *      R   : blue clique size
 *      S   : red clique size
 */

#define NT_MAX 16

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include "dSFMT.h"

#define M R-2
#define N S-2
#define URAND() dsfmt_genrand_close_open(&rstate)
typedef unsigned long ULONG;

/* Replica-specific variables ************************************************/
typedef struct
{
    int s[NV][NV];  /* edge matrix (+1=blue, -1=red) */
    int h2[NV][NV]; /* local fields */
    int energy;
} rep_t;

/* Global variables **********************************************************/
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

void init_reps(rep_t *reps, rep_t **preps)
{
    rep_t *p;
    int iT, j, k, ei, h2i;

    /* initialize each replica in a simple state with all edges blue */

    /* total energy in initial state */
    ei = binomial(NV, R);

    /* local "field" at each edge in initial state (=binomial(NV-2, R-2)) */
    h2i = R*(R-1.)/(NV*(NV-1.))*ei;

    for (iT = 0; iT < nT; iT++)
    {
        p = &reps[iT];
        preps[iT] = p;
        p->energy = ei;

        for (j = 0; j < NV; j++)
        {
            for (k = 0; k < j; k++)
            {
                p->s[j][k] = 1;
                p->h2[j][k] = h2i;
            }
        }
    }   /* end of loop over temperatures */
}

/* flip a spin and update fields */
void flip(rep_t *p, int j, int k, int delta)
{
    int l, m, isclique, a[M];

    p->energy += delta;

    /* update local field at affected edges */
    if ((p->s[j][k] *= -1) == 1)
    {
        /* loop over all M = R-2 combinations of vertices excepting j and k */
        for (i = 0; i < M; i++) a[i] = 0;
        while (a[0] < NV)
        {
            for (i=0; i < M; i++)
                printf("%d\n", a[i]);

            /* test whether this combination is a clique except for edge jk */
            isclique = true;
            for (l = 0; l < M; l++)
            {
                /* THERE'S A PROBLEM HERE */
                if (! (s[l][j] == 1 && s[l][k] == 1))
                {
                    isclique = false;
                    break;
                }
                for (m = 0; m < l; m++)
                {
                    if (! (s[l][m] == 1))
                    {
                        isclique = false;
                        break;
                    }
                }
                if (! isclique) break;
            }

            if (isclique)
            {
                for (l = 0; l < M; l++)
                {
                    for (m = 0; m < l; m++)

                }
            }

            for (i = M-1; i >= 0; i--)
            {
                while (++a[i]==j || a[i]==k);
                if (a[i] < NV) break;
                a[i] = 0;
            }
        }
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
                delta = p->s[j][k]*p->h[j][k];

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

    /* check compile-time parameters */
    assert(R > S);

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

    return EXIT_SUCCESS;
}
