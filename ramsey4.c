/*
 * To compile, run python compile2.py
 *
 * Undefined constants (computed by compile2.py)
 *      NV      : number of vertices
 *      R       : red clique size
 *      S       : blue clique size
 *      NED     : number of edges (=NV(NV-1)/2)
 *      NSGR    : number of subgraphs with R vertices (=binomial(NV, S))
 *      NSGS    : number of subgraphs with S vertices (=binomial(NV, S))
 *      NSGFER  : number of subgraphs with R vertices including a given edge
 *      NSGFES  : number of subgraphs with S vertices including a given edge
 *                  (=binomial(NV-2, S-2))
 */

#define MAX_NT          32
#define MAX_SWEEPS      10000
#define WRITE_INTERVAL  10000

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
    int h2[NED];
    int nbr[NSGR];  /* number of blue edges in each R-subgraph */
    int nbs[NSGS];  /* number of blue edges in each S-subgraph */
    int energy;     /* number of blue S-cliques and red R-cliques */
} rep_t;

/* Global variables **********************************************************/

int nedr = R*(R-1)/2;       /* number of edges in an R-subgraph */
int neds = S*(S-1)/2;       /* number of edges in an S-subgraph */
int nedsm1 = S*(S-1)/2 - 1; /* neds minus one */
int nedsm2 = S*(S-1)/2 - 2;

/* subs[ei] (subr[ei]) lists the NSGFES (NSGFER) complete S-subgraphs
 * (R-subgraphs) that include edge ei */
int *subr[NED];
int *subs[NED];
 
/* edgs[si] (edgr[si]) lists the edges of S-subgraph (R-subgraph) si */
int *edgr[NSGR];
int *edgs[NSGS];

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


void init_subgraph_table(int *sub[], int *edg[], int t,
        int nsg, int nsgfe, int nedsg)
{
    int ps[NED];
    int pe[nsg];    /* uses C99 auto dynamic arrays */
    int c[t+2];     /* array of vertices of the current subgraph */
    int ei, si;     /* edge index, subgraph index */
    int j, k;

    for (j = 0; j < NED; j++)
    {
        sub[j] = (int*) malloc(nsgfe * sizeof(int));
        ps[j] = 0;
    }

    for (j = 0; j < nsg; j++)
    {
        edg[j] = (int*) malloc(nedsg * sizeof(int));
        pe[j] = 0;
    }

    /* 
     * iterate over all subgraphs with t vertices
     */

    /*
     * algorithm to generate combinations adapted from Algorithm L in Knuth's
     * Art of Computer Programming Vol. 4, Fasc. 3 (all-caps labels
     * correspond to labels in the book)
     */

    /* INITIALIZE */
    si = 0;
    c[t] = NV;
    c[t+1] = 0;
    for (j = 0; j < t; j++) c[j] = j;

    while (1)
    {
        /*
         * VISIT combination c_1 c_2 ... c_t
         * (algorithm guarantees that c_1 < c_2 < ... < c_t)
         */

        /* iterate over edges in this subgraph */
        for (k = 0; k < t; k++)
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
        if (j == t) break;

        c[j]++;
    }
}

void free_subgraph_table(int *sub[], int *edg[], int nsg)
{
    int i;

    for (i = 0; i < NED; i++)
        free(sub[i]);
    for (i = 0; i < nsg; i++)
        free(edg[i]);
}

void update(int ei, int sp[], int nbr[], int nbs[], int h2[])
{
    int si, j, ej, nbf;

    if (sp[ei] == 1)
    {
        for (si = 0; si < NSGFES; si++)
        {
            nbf = nbs[subs[ei][si]] += 1;

            if (nbf == neds)    /* created blue clique */
            {
                h2[ei] += 1;
                for (j = 0; j < neds; j++) h2[edgs[subs[ei][si]][j]] -= 1;
            }
            else if (nbf == nedsm1) /* created incomplete blue clique */
            {
                for (j = 0; j < neds; j++)
                {
                    ej = edgs[subs[ei][si]][j];
                    if (sp[ej] == -1) { h2[ej] -= 1; break; }
                }
            }
        }

        for (si = 0; si < NSGFER; si++)
        {
            nbf = nbr[subr[ei][si]] += 1;

            if (nbf == 1)   /* destroyed red clique */
            {
                h2[ei] += 1;
                for (j = 0; j < nedr; j++) h2[edgr[subr[ei][si]][j]] -= 1;
            }
            else if (nbf == 2)   /* destroyed incomplete red clique */
            {
                for (j = 0; j < nedr; j++)
                {
                    ej = edgr[subr[ei][si]][j];
                    if (sp[ej] == 1 && ej != ei) { h2[ej] -= 1; break; }
                }
            }
        }
    }
    else
    {
        for (si = 0; si < NSGFER; si++)
        {
            nbf = nbr[subr[ei][si]] -= 1;

            if (nbf == 0)   /* created a red clique */
            {
                h2[ei] -= 1;
                for (j = 0; j < nedr; j++) h2[edgr[subr[ei][si]][j]] += 1;
            }
            else if (nbf == 1)  /* created an incomplete red clique */
            {
                for (j = 0; j < nedr; j++)
                {
                    ej = edgr[subr[ei][si]][j];
                    if (sp[ej] == 1) { h2[ej] += 1; break; }
                }
            }
        }

        for (si = 0; si < NSGFES; si++)
        {
            nbf = nbs[subs[ei][si]] -= 1;

            if (nbf == nedsm1)  /* destroyed a blue clique */
            {
                h2[ei] -= 1;
                for (j = 0; j < neds; j++) h2[edgs[subs[ei][si]][j]] += 1;
            }
            else if (nbf == nedsm2) /* destroyed an incomplete blue clique */
            {
                for (j = 0; j < neds; j++)
                {
                    ej = edgs[subs[ei][si]][j];
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
        p->energy = NSGS;
        nswaps[it] = 0;

        for (j = 0; j < NSGR; j++) p->nbr[j] = nedr;
        for (j = 0; j < NSGS; j++) p->nbs[j] = neds;

        for (j = 0; j < NED; j++)
        {
            p->sp[j] = 1;
            p->h2[j] = -NSGFES;
        }

        /* randomize spins */
        for (j = 0; j < NED; j++)
        {
            if (URAND() < 0.5)
            {
                p->sp[j] = -1;
                p->energy += p->h2[j];
                update(j, p->sp, p->nbr, p->nbs, p->h2);
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
                update(j, p->sp, p->nbr, p->nbs, p->h2);
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
    fprintf(fp, "%d\n", R);
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
#ifndef NOTIME
    start = clock();
#endif

    while (! done && nsweeps < MAX_SWEEPS)
    {
        sweep();
        nsweeps++;

        temper();

        if (nsweeps % WRITE_INTERVAL == 0)
        {
            sprintf(filename, "%d-%d-%d_%d.bin",
                   R, S, NV, rseed);
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
                sprintf(filename, "%d-%d-%d_%d.graph", R, S, NV, rseed);
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

    init_subgraph_table(subr, edgr, R, NSGR, NSGFER, nedr);
    init_subgraph_table(subs, edgs, S, NSGS, NSGFES, neds);
    init_replicas();

    if (argc == 4) load_state(argv[3]);

    run();
    
    free_subgraph_table(subr, edgr, NSGR);
    free_subgraph_table(subs, edgs, NSGS);

    return EXIT_SUCCESS;
}
