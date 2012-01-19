/*
 * To compile, run python compile2.py
 *
 * Undefined constants (computed by compile2.py)
 *      NV      : number of vertices
 *      R       : red clique size
 *      S       : blue clique size
 *      NSGR    : number of subgraphs with R vertices (=binomial(NV, S))
 *      NSGS    : number of subgraphs with S vertices (=binomial(NV, S))
 *      NSGFER  : number of subgraphs with R vertices including a given edge
 *      NSGFES  : number of subgraphs with S vertices including a given edge
 *                  (=binomial(NV-2, S-2))
 */

#define MAX_NT          32
#define MAX_SWEEPS      100000
#define WRITE_INTERVAL  100000

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "dSFMT.h"

#define URAND() dsfmt_genrand_close_open(&rstate)
typedef unsigned long ULONG;

/* Replica-specific variables ************************************************/
typedef struct
{
    int sp[NV][NV];  
    /* for j < k, sp[j][k] = 1 if edge (j, k) is blue, -1 if edge is red */

    int nbr[NSGR];  /* number of blue edges in each R-subgraph */
    int nbs[NSGS];  /* number of blue edges in each S-subgraph */
    int energy;     /* number of blue S-cliques and red R-cliques */
} rep_t;

/* Global variables **********************************************************/

int nedr;    /* number of edges in an R-subgraph (=S(S-1)/2) */
int neds;    /* number of edges in an S-subgraph (=S(S-1)/2) */
int nedrm1;  /* nedr minus one */
int nedsm1;  /* neds minus one */

int *subr[NV][NV];
int *subs[NV][NV];
/*
 * for j < k, subs[j][k] (subr[j][k]) is an array of length NSGFE containing
 * the labels of all subgraphs with S (R) vertices that include the edge (j, k)
 */

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


void init_subgraph_table(int *sub[NV][NV], int t, int nsgfe)
{
    int len[NV][NV];    /* positions in subgraph arrays */
    int c[t+2];         /* array of vertices of the current subgraph */
    int id;             /* subgraph label */
    int j, k;
    int cj, ck;

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
     * iterate over all subgraphs with t vertices
     */

    /*
     * algorithm to generate combinations adapted from Algorithm L in Knuth's
     * Art of Computer Programming Vol. 4, Fasc. 3 (all-caps labels
     * correspond to labels in the book)
     */

    /* INITIALIZE */
    id = 0;
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
        if (j == t) break;

        c[j]++;
    }
}

void free_subgraph_table(int *sub[NV][NV])
{
    int j, k;

    for (k = 0; k < NV; k++)
        for (j = 0; j < k; j++)
            free(sub[j][k]);
}

int flip_energy(int sp, int *subr, int *subs, int *nbr, int *nbs)
{
    int i, delta;

    delta = 0;

    if (sp == 1)
    {
        for (i = 0; i < NSGFER; i++)
            if (nbr[subr[i]] == 1) delta++;
        for (i = 0; i < NSGFES; i++)
            if (nbs[subs[i]] == neds) delta--;
    }
    else
    {
        for (i = 0; i < NSGFER; i++)
            if (nbr[subr[i]] == 0) delta--;
        for (i = 0; i < NSGFES; i++)
            if (nbs[subs[i]] == nedsm1) delta++;
    }

    return delta;
}

void update_nb(int sp, int *subr, int *subs, int *nbr, int *nbs)
{
    int i;

    for (i = 0; i < NSGFER; i++) nbr[subr[i]] += sp;
    for (i = 0; i < NSGFES; i++) nbs[subs[i]] += sp;
}

/* initialize each replica in a random configuration */
void init_replicas()
{
    rep_t *p;
    int it, j, k;

    for (it = 0; it < nt; it++)
    {
        ri[it] = it;
        p = &reps[it];
        p->energy = NSGS;
        nswaps[it] = 0;

        for (j = 0; j < NSGR; j++) p->nbr[j] = nedr;
        for (j = 0; j < NSGS; j++) p->nbs[j] = neds;

        for (k = 0; k < NV; k++)
        {
            for (j = 0; j < k; j++)
            {
                if (URAND() > 0.5) p->sp[j][k] = 1;
                else
                {
                    p->sp[j][k] = -1;
                    p->energy += flip_energy(1, subr[j][k], subs[j][k], p->nbr, p->nbs);
                    update_nb(-1, subr[j][k], subs[j][k], p->nbr, p->nbs);
                }
            }
        }
    }
}

void sweep()
{
    rep_t *p;
    int it, j, k, delta;

    for (it = 0; it < nt; it++)
    {
        p = &reps[ri[it]];
        for (k = 0; k < NV; k++)
        {
            for (j = 0; j < k; j++)
            {
                /* compute energy difference of flip */
                delta = flip_energy(p->sp[j][k], subr[j][k], subs[j][k], p->nbr, p->nbs); 

                /* flip with Metropolis probability */
                if (delta <= 0 || URAND() < exp(mbeta[it]*delta))
                {
                    p->energy += delta;
                    p->sp[j][k] *= -1;
                    update_nb(p->sp[j][k], subr[j][k], subs[j][k], p->nbr, p->nbs);
                }
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

void save_graph(int sp[NV][NV], char filename[])
{
    FILE *fp;
    int j, k;

    fp = fopen(filename, "w");
    fprintf(fp, "%d\n", NV);
    fprintf(fp, "%d\n", R);
    fprintf(fp, "%d\n", S);
    
    for (k = 0; k < NV; k++)
        for (j = 0; j < k; j++)
            fprintf(fp, "%d\n", (sp[j][k] == 1) ? 1 : 0);

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

    printf("%3d | ", min);
    for (it = 0; it < nt; it++)
        printf("%3d ", reps[ri[it]].energy);
    for (it = 1; it < nt; it++)
        printf("%3.2f ", (float) nswaps[it]/nsweeps);
    printf("\n");
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

        /*if (nsweeps % 100 == 0) print_status();*/

        if (nsweeps % WRITE_INTERVAL == 0)
        {
            sprintf(filename, "%d-%d-%d_%d.bin",
                   R, S, NV, nsweeps/WRITE_INTERVAL);
            save_state(filename);
            printf("State saved to %s\n", filename);
        }

        for (it = 0; it < nt; it++)
        {
            p = &reps[ri[it]];
            if (p->energy < min)
            {
                min = p->energy;
                print_status();

                if (p->energy == 0)
                {
                    sprintf(filename, "%d-%d-%d_0.graph", R, S, NV);
                    save_graph(reps[ri[it]].sp, filename);

                    printf("Found zero-energy ground state!\n");
                    printf("Graph saved to %s\n", filename);
                    printf("N_sweeps = %d\n", nsweeps);

                    done = 1;
                    break;
                }
            }
        }
    }

    sprintf(filename, "%d-%d-%d_%d.graph", R, S, NV, reps[0].energy);
    save_graph(reps[ri[it]].sp, filename);
    printf("Graph saved to %s\n", filename);
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

    nedr = R*(R-1)/2;
    neds = S*(S-1)/2;
    nedrm1 = nedr - 1;
    nedsm1 = neds - 1;

    init_subgraph_table(subr, R, NSGFER);
    init_subgraph_table(subs, S, NSGFES);
    init_replicas();

    if (argc == 4) load_state(argv[3]);

    run();
    
    free_subgraph_table(subr);
    free_subgraph_table(subs);

    return EXIT_SUCCESS;
}
