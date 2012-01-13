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

#define URAND() dsfmt_genrand_close_open(&rstate)
typedef unsigned long ULONG;

/* Replica-specific variables ************************************************/
typedef struct
{
    int s[NV][NV];  /* edge (spin) matrix (+1=blue, -1=red) */
    int h2[NV][NV]; /* field */
    int *m;
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
    int len[NV][NV];
    int c[S+2];
    int i, j, ci, cj, id;

    nsg = binomial(NV, S);
    nsg_fe = S*(S-1.)/(NV*(NV-1.))*nsg; /* (=binomial(NV-2, S-2)) */
    /* debug */
    /* assert(nsg_fe==binomial(NV-2,S-2)); */

    /* initialize subgraph lists */
    for (i = 0; i < NV; i++)
    {
        for (j = 0; j < i; j++)
        {
            sub[i][j] = (int*) malloc(nsg_fe * sizeof(int));
            len[i][j] = 0;
        }
    }

    /* 
     * iterate over all subgraphs with S vertices
     * ------------------------------------------------------------------------
     * algorithm to generate combinations adapted from Algorithm L in Knuth's
     * Art of Computer Programming Vol. 4, Fasc. 3 (all-caps labels
     * correspond to labels in the book)
     */

    /* INITIALIZE */
    c[S] = NV;
    c[S+1] = 0;
    for (i = 0; i < S; i++) c[i] = i;

    id = 0;

    while (1)
    {
        /* VISIT combination c_1 c_2 ... c_S */
        /* (algorithm guarantees that c_1 < c_2 < ... < c_S) */
        /* iterate over edges in this subgraph */
        for (i = 0; i < S; i++)
        {
            for (j = 0; j < i; j++)
            {
                ci = c[i];
                cj = c[j];

                /* add subgraph to list for edge (i, j) */
                sub[ci][cj][len[ci][cj]++] = id;
            }
        }

        /* FIND j */
        for (j = 0; j <= S; j++)
        {
            if (c[j] + 1 == c[j+1]) c[j] = j;
            else break;
        }

        /* DONE? */
        if (j == S) break;

        c[j]++;
        id++;   /* increment label for next combination */
    }
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
        p->m = (int*) malloc(nsg * sizeof(int));

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

void free_reps()
{
    int iT;

    for (iT = 0; iT < nT; iT++)
        free(reps[iT].m);
}

/* flip a spin and update fields */
void flip(rep_t *p, int j, int k, int delta)
{
    p->energy += delta;

    /* update local field at affected edges */

    if ((p->s[j][k] *= -1) == 1)
    {
        for (i = 0; i < nsg; i++)
            if (++p->m[sub[j][k][i]] == S)
            {
                /* flip has created a new blue clique */
            }
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
