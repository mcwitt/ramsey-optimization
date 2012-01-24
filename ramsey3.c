/*
 * To compile, run python compile3.py
 *
 * Undefined constants (computed by compile.py)
 *      NV      : number of vertices
 *      S       : clique size
 *      NED     : number of edges (=NV(NV-1)/2)
 *      NSG     : number of subgraphs with S vertices (=binomial(NV, S))
 *      NSGFE   : number of subgraphs with S vertices including a given edge
 *                  (=binomial(NV-2, S-2))
 */

#define MAX_NT          32
#define MAX_SWEEPS      100
#define WRITE_INTERVAL  1000

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "dSFMT.h"

#define URAND() dsfmt_genrand_close_open(&rstate)

/* Replica-specific variables ************************************************/
typedef struct
{
    int sp[NED];
    int h2[NED];    /* local field at each edge (doubled) */
    int nb[NSG];    /* number of blue edges in each S-subgraph */
    int energy;     /* number of blue S-cliques and red S-cliques */
} rep_t;

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

void init_subgraph_table()
{
    int ps[NED];
    int pe[NSG];
    int c[S+2];     /* array of vertices of the current subgraph */
    int ei, si;     /* edge index, subgraph index */
    int j, k;

    for (j = 0; j < NED; j++)
    {
        sub[j] = (int*) malloc(NSGFE * sizeof(int));
        ps[j] = 0;
    }

    for (j = 0; j < NSG; j++)
    {
        edg[j] = (int*) malloc(neds * sizeof(int));
        pe[j] = 0;
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
    si = 0;
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
                ei = c[k]*(c[k]-1)/2 + c[j];

                /*
                 * add subgraph si to list for edge ei,
                 * and edge ei to list for subgraph si
                 */

                sub[ei][ps[ei]++] = si;
                edg[si][pe[si]++] = ei;
            }
        }

        si++;   /* finished with this subgraph, increment label */

        /* FIND j */
        j = 0;
        while (c[j] + 1 == c[j+1]) { c[j] = j; j++; }

        /* DONE? */
        if (j == S) break;

        c[j]++;
    }
}

void free_subgraph_table()
{
    int i;

    for (i = 0; i < NED; i++)
        free(sub[i]);
    for (i = 0; i < NSG; i++)
        free(edg[i]);
}

void update(int ei, int sp[], int nb[], int h2[])
{
    int si, j, ej, nbf;

    if (sp[ei] == 1)
    {
        for (si = 0; si < NSGFE; si++)
        {
            nbf = nb[sub[ei][si]] += 1;
            if (nbf == neds || nbf == 1)
            {
                /* completed a blue clique or destroyed a red clique */
                h2[ei] += 1;
                for (j = 0; j < neds; j++)
                    h2[edg[sub[ei][si]][j]] -= 1;
            }
            else if (nbf == nedsm1)
            {
                /* completed blue clique except for one red edge--
                 * update field at red edge */
                for (j = 0; j < neds; j++)
                {
                    ej = edg[sub[ei][si]][j];
                    if (sp[ej] == -1) { h2[ej] -= 1; break; }
                }
            }
            else if (nbf == 2)
            {
                /* destroyed an almost-complete red clique--
                 * update field at the preexisting blue edge */
                for (j = 0; j < neds; j++)
                {
                    ej = edg[sub[ei][si]][j];
                    if (sp[ej] == 1 && ej != ei) { h2[ej] -= 1; break; }
                }
            }
        }
    }
    else    /* analagous procedure if the new edge is red... */
    {
        for (si = 0; si < NSGFE; si++)
        {
            nbf = nb[sub[ei][si]] -= 1;
            if (nbf == 0 || nbf == nedsm1)
            {
                h2[ei] -= 1;
                for (j = 0; j < neds; j++)
                    h2[edg[sub[ei][si]][j]] += 1;
            }
            else if (nbf == 1)
            {
                for (j = 0; j < neds; j++)
                {
                    ej = edg[sub[ei][si]][j];
                    if (sp[ej] == 1) { h2[ej] += 1; break; }
                }
            }
            else if (nbf == nedsm2)
            {
                for (j = 0; j < neds; j++)
                {
                    ej = edg[sub[ei][si]][j];
                    if (sp[ej] == -1 && ej != ei) { h2[ej] += 1; break; }
                }
            }
        }
    }
}

#ifdef DEBUG
#include "debug.c"
#endif

/* initialize each replica in a random configuration */
void init_replicas()
{
    rep_t *p;
    int it, j;

    for (it = 0; it < nt; it++)
    {
        ri[it] = it;
        p = &reps[it];
        p->energy = NSG;
        nswaps[it] = 0;

        for (j = 0; j < NSG; j++) p->nb[j] = neds;

        for (j = 0; j < NED; j++)
        {
            p->sp[j] = 1;
            p->h2[j] = -NSGFE;
        }

        /* randomize spins */
        for (j = 0; j < NED; j++)
        {
            if (URAND() < 0.5)
            {
                p->sp[j] = -1;
                p->energy += p->h2[j];
                update(j, p->sp, p->nb, p->h2);
            }
        }
    }
}

void sweep()
{
    rep_t *p;
    int it, j, delta;

    for (it = 0; it < nt; it++)
    {
        p = &reps[ri[it]];

        for (j = 0; j < NED; j++)
        {
            /* compute energy difference of flip */
            delta = p->h2[j]*p->sp[j];

            /* flip with Metropolis probability */
            if (delta <= 0 || URAND() < exp(mbeta[it]*delta))
            {
                p->sp[j] *= -1;
                p->energy += delta;
                update(j, p->sp, p->nb, p->h2);
            }
        }
#ifdef DEBUG
        assert(debug_energy(p->sp) == p->energy);
#endif
    }   /* end of loop over temperatures */
}

void temper()
{
    double logar;
    int it, copy;

    for (it = 1; it < nt; it++)
    {
        logar = (reps[ri[it-1]].energy - reps[ri[it]].energy)
            * (mbeta[it] - mbeta[it-1]);

        if (URAND() < exp(logar))
        {
            /* do PT swap */
            copy = ri[it-1];
            ri[it-1] = ri[it];
            ri[it] = copy;
            nswaps[it]++;
        }
    }
}

void save_graph(int sp[NED], char filename[])
{
    FILE *fp;
    int i;

    fp = fopen(filename, "w");
    fprintf(fp, "%d\n", NV);
    fprintf(fp, "%d\n", S);
    fprintf(fp, "%d\n", S);
    
    for (i = 0; i < NED; i++)
        fprintf(fp, "%d\n", (sp[i] == 1) ? 1 : 0);

    fclose(fp);
}

void save_state(char filename[])
{
    FILE *fp;

    fp = fopen(filename, "w");
    fwrite(ri, sizeof(int), nt, fp);
    fwrite(reps, sizeof(rep_t), nt, fp);
    fclose(fp);
}

void load_state(char filename[])
{
    FILE *fp;

    fp = fopen(filename, "r");
    assert(fread(ri, sizeof(int), nt, fp) == nt);
    assert(fread(reps, sizeof(rep_t), nt, fp) == nt);
    fclose(fp);
}

void print_status()
{
    int it;
#ifndef NOTIME
    int trun;
#endif

    printf("\n");
    printf("min. energy     : %d\n", min);
    printf("# of sweeps     : %d\n", nsweeps);
#ifndef NOTIME
    trun = (clock() - start)/CLOCKS_PER_SEC;
    printf("time running    : %d seconds\n", trun);
    printf("sweep rate      : %.2f / s\n", (float) nsweeps/trun);
#endif
    for (it = 0; it < nt; it++)
        printf("%5d ", reps[ri[it]].energy);
    printf("\n");
    for (it = 0; it < nt; it++)
        printf("%5.2f ", T[it]);
    printf("\n   ");
    for (it = 1; it < nt; it++)
        printf("%5.2f ", (float) nswaps[it]/nsweeps);
    printf("\n");
    fflush(stdout);
}

void run()
{
    rep_t *p;
    char filename[256];
    int it, done;

    min = INT_MAX;
    nsweeps = 0;
    done = 0;

    while (! done && nsweeps < MAX_SWEEPS)
    {
        sweep();
        nsweeps++;

        temper();

        if (nsweeps % WRITE_INTERVAL == 0)
        {
            sprintf(filename, "%d-%d-%d_%d.bin",
                    S, S, NV, nsweeps/WRITE_INTERVAL);
            save_state(filename);
            printf("state saved to %s\n", filename);
        }

        for (it = 0; it < nt; it++)
        {
            p = &reps[ri[it]];
            if (p->energy < min)
            {
                min = p->energy;
                print_status();
                sprintf(filename, "%d-%d-%d_%d.graph", S, S, NV, rseed);
                save_graph(p->sp, filename);

                if (p->energy == 0) { done = 1; break; }
            }
        }
    }

    print_status();
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

    init_subgraph_table();
    init_replicas();

    if (argc == 4) load_state(argv[3]);

    run();
    
    free_subgraph_table();

    return EXIT_SUCCESS;
}
