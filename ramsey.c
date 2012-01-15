/*
 * Undefined constants
 *      NV  : number of vertices
 *      S   : clique size
 */

#define NT_MAX 32

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
    int *nb;
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
    int j, k, cj, ck, id;

    nsg = binomial(NV, S);
    nsg_fe = S*(S-1.)/(NV*(NV-1.))*nsg; /* (=binomial(NV-2, S-2)) */
    /* debug */
    /* assert(nsg_fe==binomial(NV-2,S-2)); */

    /* initialize subgraph lists */
    for (k = 0; k < NV; k++)
    {
        for (j = 0; j < k; j++)
        {
            sub[j][k] = (int*) malloc(nsg_fe * sizeof(int));
            len[j][k] = 0;
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
    for (j = 0; j < S; j++) c[j] = j;

    id = 0;

    while (1)
    {
        /* VISIT combination c_1 c_2 ... c_S */
        /* (algorithm guarantees that c_1 < c_2 < ... < c_S) */
        /* iterate over edges in this subgraph */
        for (k = 0; k < S; k++)
        {
            for (j = 0; j < k; j++)
            {
                cj = c[j];
                ck = c[k];

                /* add subgraph to list for edge (i, j) */
                sub[cj][ck][len[cj][ck]++] = id;
            }
        }

        id++;   /* increment label for next combination */

        /* FIND j */
        for (j = 0; j <= S; j++)
        {
            if (c[j] + 1 == c[j+1]) c[j] = j;
            else break;
        }

        /* DONE? */
        if (j == S) break;

        c[j]++;
    }
}

void free_globals()
{
    int j, k;

    for (k = 0; k < NV; k++)
        for (j = 0; j < k; j++)
            free(sub[j][k]);
}

/* initialize each replica with all spins +1 (i.e. all edges blue) */
void flip(rep_t *p, int j, int k, int delta);
int flip_energy(int s, int *sub, int *nb);
void init_reps()
{
    rep_t *p;
    int iT, j, k;

    for (iT = 0; iT < nT; iT++)
    {
        p = &reps[iT];
        preps[iT] = p;
        p->energy = nsg;
        p->nb = (int*) malloc(nsg * sizeof(int));

        for (j = 0; j < nsg; j++)
            p->nb[j] = S;

        for (k = 0; k < NV; k++)
            for (j = 0; j < k; j++)
            {
                p->s[j][k] = 1;
                if (URAND() < 0.5)
                    flip(p, j, k, flip_energy(1, sub[j][k], p->nb));
            }

    }
}

void free_reps()
{
    int iT;

    for (iT = 0; iT < nT; iT++)
        free(reps[iT].nb);
}

int flip_energy(int s, int *sub, int *nb)
{
    int i, nbi, delta;

    delta = 0;

    if (s == 1)
    {
        for (i = 0; i < nsg_fe; i++)
        {
            nbi = nb[sub[i]];
            if (nbi == S) delta--;
            else if (nbi == 1) delta++;
        }
    }
    else
    {
        for (i = 0; i < nsg_fe; i++)
        {
            nbi = nb[sub[i]];
            if (nbi == 0) delta--;
            else if (nbi == S - 1) delta++;
        }
    }

    return delta;
}

void flip(rep_t *p, int j, int k, int delta)
{
    int i;

    p->energy += delta;
    p->s[j][k] *= -1;

    for (i = 0; i < nsg_fe; i++)
        p->nb[sub[j][k][i]] += p->s[j][k];
}

void sweep()
{
    rep_t *p;
    int iT, j, k, delta;

    for (iT = 0; iT < nT; iT++)
    {
        p = preps[iT];
        for (k = 0; k < NV; k++)
        {
            for (j = 0; j < k; j++)
            {
                /* compute energy difference of flip */
                delta = flip_energy(p->s[j][k], sub[j][k], p->nb);

                /* flip with Metropolis probability */
                if (delta <= 0 || URAND() < exp(mbeta[iT]*delta))
                    flip(p, j, k, delta);
            }
        }
    }   /* end of loop over temperatures */
}

void temper()
{
    rep_t *pl, *pr;
    double logar;
    int iT;

    for (iT = 1; iT < nT; iT++)
    {
        pl = preps[iT-1];
        pr = preps[iT];

        /* compute log of acceptance ratio for swap */
        logar = (mbeta[iT]-mbeta[iT-1])*(pl->energy-pr->energy);

        if (URAND() < exp(logar))
        {
            /* do PT swap (swap entries iT and iT-1 in preps) */
            preps[iT-1] = pr;
            preps[iT] = pl;
            nswaps[iT]++;
        }
    }
}

void run()
{
    int iT, nsweeps, done;

    nsweeps = 0;
    done = 0;

    while (! done)
    {
        sweep();

        if (nsweeps % 10 == 0)
        {
            temper();

            if (nsweeps % 20 == 0)
            {
            for (iT=0; iT<nT; iT++)
                printf("%3d ", preps[iT]->energy);
            for (iT=1; iT<nT; iT++)
                printf("%3.2f ", (float) nswaps[iT]/nsweeps);

            printf("\n");
            }

        }


        for (iT = 0; iT < nT; iT++)
        {
            if (preps[iT]->energy == 0)
            {
                printf("Found zero-energy ground state\n");
                printf("%d\n", nsweeps);
                done = 1;
                break;
            }
        }

        nsweeps++;
    }
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


    /* init random number generator */
    rseed = atoi(argv[2]);
    dsfmt_init_gen_rand(&rstate, rseed);

    /* read temperatures from input file */
    nT = 0;
    infile = fopen(argv[1], "r");
    while (fscanf(infile, "%lf", &t) != EOF)
    {
        T[nT] = t;
        mbeta[nT] = -1./t;
        nT++;
    }
    
    assert(1 < nT && nT < NT_MAX);

    init_globals();
    init_reps();

    run();
    
    free_reps();
    free_globals();

    return EXIT_SUCCESS;
}
