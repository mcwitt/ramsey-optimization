/*
 * Undefined constants
 *      NV  : number of vertices
 *      S   : clique size
 */

#define NT_MAX 24

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "dSFMT.h"

#define URAND() dsfmt_genrand_close_open(&rstate)
typedef unsigned long ULONG;

/* Replica-specific variables ************************************************/
typedef struct
{
    int s[NV][NV];  
    /* for j < k, s[j][k] = 1 if edge (j, k) is blue, -1 if edge is red */

    int *nb;        /* number of blue edges in each S-subgraph */
    int energy;     /* number of blue S-cliques and red S-cliques */
} rep_t;

/* Global variables **********************************************************/

int ned;    /* number of edges in an S-subgraph (=S(S-1)/2) */
int nedm1;  /* ned minus one */
int nsg;    /* number of subgraphs with S vertices (=binomial(N, S)) */
int nsgfe;  /* number of subgraphs including a given edge */

int *sub[NV][NV];   
/*
 * for j < k, sub[j][k] is an array of length nsgfe containing the labels of
 * all subgraphs with S vertices that include the edge (j, k)
 */

rep_t reps[NT_MAX]; /* storage for parallel tempering (PT) replicas */

/* pointers to PT replicas in order of increasing temperature */
rep_t *preps[NT_MAX];   

int nT;                 /* number of PT copies */
int nswaps[NT_MAX];     /* number of swaps between each pair of temperatures */
double T[NT_MAX];       /* array of temperatures */
double mbeta[NT_MAX];   /* negative inverse temperatures */

dsfmt_t rstate; /* state of random number generator (RNG) */
uint32_t rseed; /* seed used to initialize RNG */

/* include debugging functions */
#include "debug.c"

/* binomial coefficient */
ULONG binomial(ULONG n, ULONG k)
{
    ULONG r = 1, d = n - k; 

    /*if (k > n) return 0;*/

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
    int len[NV][NV];    /* positions in subgraph arrays */
    int c[S+2];         /* array of vertices of the current subgraph */
    int id;             /* subgraph label */
    int j, k;
    int cj, ck;

    ned = S*(S-1)/2;
    nedm1 = ned - 1;
    nsg = binomial(NV, S);
    nsgfe = S*(S-1.)/(NV*(NV-1.))*nsg; /* (=binomial(NV-2, S-2)) */

    /* initialize subgraph arrays */
    for (k = 0; k < NV; k++)
    {
        for (j = 0; j < k; j++)
        {
            sub[j][k] = (int*) malloc(nsgfe * sizeof(int));
            len[j][k] = 0;
        }
    }

    /* 
     * iterate over all subgraphs with S vertices
     */

    /*
     * algorithm to generate combinations adapted from Algorithm L in Knuth's
     * Art of Computer Programming Vol. 4, Fasc. 3 (all-caps labels
     * correspond to labels in the book)
     */

    /* INITIALIZE */
    id = 0;
    c[S] = NV;
    c[S+1] = 0;
    for (j = 0; j < S; j++) c[j] = j;

    while (1)
    {
        /*
         * VISIT combination c_1 c_2 ... c_S
         * (algorithm guarantees that c_1 < c_2 < ... < c_S)
         */

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

        id++;   /* finished with this subgraph, increment label */

        /* FIND j */
        j = 0;
        while (c[j] + 1 == c[j+1]) { c[j] = j; j++; }

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

int flip_energy(int s, int *sub, int *nb)
{
    int i, nbi, delta;

    delta = 0;

    if (s == 1)
    {
        for (i = 0; i < nsgfe; i++)
        {
            nbi = nb[sub[i]];
            if (nbi == ned) delta--;
            else if (nbi == 1) delta++;
        }
    }
    else
    {
        for (i = 0; i < nsgfe; i++)
        {
            nbi = nb[sub[i]];
            if (nbi == 0) delta--;
            else if (nbi == nedm1) delta++;
        }
    }

    return delta;
}

void update_nb(int s, int *sub, int *nb)
{
    int i;

    for (i = 0; i < nsgfe; i++) nb[sub[i]] += s;
}

/* initialize each replica with a random configuration */
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

        for (j = 0; j < nsg; j++) p->nb[j] = ned;

        for (k = 0; k < NV; k++)
        {
            for (j = 0; j < k; j++)
            {
                if (URAND() > 0.5) p->s[j][k] = 1;
                else
                {
                    p->s[j][k] = -1;
                    p->energy += flip_energy(1, sub[j][k], p->nb);
                    update_nb(-1, sub[j][k], p->nb);
                }
            }
        }
    }
}

void free_reps()
{
    int iT;

    for (iT = 0; iT < nT; iT++)
        free(reps[iT].nb);
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
                {
                    p->energy += delta;
                    update_nb(p->s[j][k] *= -1, sub[j][k], p->nb);
                }
            }
        }
        assert(debug_energy(p->s) == p->energy);
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

void write(int s[NV][NV], char filename[])
{
    FILE *fp;
    int j, k;

    fp = fopen(filename, "w");

    fprintf(fp, "%d\n", NV);
    fprintf(fp, "%d\n", S);
    fprintf(fp, "%d\n", S);

    for (k = 0; k < NV; k++)
        for (j = 0; j < k; j++)
            fprintf(fp, "%d\n", (s[j][k] == 1) ? 1 : 0);

    fclose(fp);
}

void run()
{
    int iT, nsweeps, done;

    nsweeps = 0;
    done = 0;

    while (! done)
    {
        sweep();
        temper();

        if (nsweeps % 10 == 0)
        {

            for (iT=0; iT<nT; iT++)
                printf("%3d ", preps[iT]->energy);
            for (iT=1; iT<nT; iT++)
                printf("%3.2f ", (float) nswaps[iT]/nsweeps);

            printf("\n");
        }


        for (iT = 0; iT < nT; iT++)
        {
            if (preps[iT]->energy == 0)
            {
                printf("Found zero-energy ground state\n");
                printf("%d\n", nsweeps);
                write(preps[iT]->s, "zero.data");
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
