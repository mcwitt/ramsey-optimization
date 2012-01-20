/*
 * To compile, run python compile.py
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
#define MAX_SWEEPS      1000
#define WRITE_INTERVAL  1000

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "dSFMT.h"

#define URAND() dsfmt_genrand_close_open(&rstate)

/* computes index of edge (j, k), j < k */
#define EDGE(j, k) (k*(k-1)/2 + j)

typedef unsigned long ULONG;

/* Replica-specific variables ************************************************/
typedef struct
{
    int sp[NED];

    int nb[NSG];    /* number of blue edges in each S-subgraph */
    int energy;     /* number of blue S-cliques and red S-cliques */
} rep_t;

/* Global variables **********************************************************/

int ned;    /* number of edges in an S-subgraph (=S(S-1)/2) */
int nedm1;  /* ned minus one */

int *sub[NED];
/*
 * for j < k, sub[j][k] is an array of length NSGFE containing the labels of
 * all subgraphs with S vertices that include the edge (j, k)
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


void init_subgraph_table()
{
    int len[NED];       /* positions in subgraph arrays */
    int c[S+2];         /* array of vertices of the current subgraph */
    int ei, si;         /* edge index, subgraph index */
    int j, k;

    /* initialize subgraph arrays */
    for (j = 0; j < NED; j++)
    {
        sub[j] = (int*) malloc(NSGFE * sizeof(int));
        len[j] = 0;
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
                ei = EDGE(c[j], c[k]);

                /* add subgraph to list for edge (i, j) */
                sub[ei][len[ei]++] = si;
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
}

int flip_energy(int sp, int *sub, int *nb)
{
    int i, nbi, delta;

    delta = 0;

    if (sp == 1)
    {
        for (i = 0; i < NSGFE; i++)
        {
            nbi = nb[sub[i]];
            if (nbi == ned) delta--;
            else if (nbi == 1) delta++;
        }
    }
    else
    {
        for (i = 0; i < NSGFE; i++)
        {
            nbi = nb[sub[i]];
            if (nbi == 0) delta--;
            else if (nbi == nedm1) delta++;
        }
    }

    return delta;
}

void update_nb(int sp, int *sub, int *nb)
{
    int i;

    for (i = 0; i < NSGFE; i++) nb[sub[i]] += sp;
}

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

        for (j = 0; j < NSG; j++) p->nb[j] = ned;

        for (j = 0; j < NED; j++)
        {
            if (URAND() > 0.5) p->sp[j] = 1;
            else
            {
                p->sp[j] = -1;
                p->energy += flip_energy(1, sub[j], p->nb);
                update_nb(-1, sub[j], p->nb);
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
            delta = flip_energy(p->sp[j], sub[j], p->nb);

            /* flip with Metropolis probability */
            if (delta <= 0 || URAND() < exp(mbeta[it]*delta))
            {
                p->energy += delta;
                update_nb(p->sp[j] *= -1, sub[j], p->nb);
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

    printf("\n");
    printf("E_min    = %d\n", min);
    printf("N_sweeps = %d\n", nsweeps);
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

        /*if (nsweeps % 100 == 0) print_status();*/

        if (nsweeps % WRITE_INTERVAL == 0)
        {
            sprintf(filename, "%d-%d-%d_%d.bin",
                    S, S, NV, nsweeps/WRITE_INTERVAL);
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
                sprintf(filename, "%d-%d-%d_min.graph", S, S, NV);
                save_graph(p->sp, filename);

                if (p->energy == 0)
                {
                    printf("Found zero-energy ground state!\n");
                    printf("N_sweeps = %d\n", nsweeps);
                    done = 1;
                    break;
                }
            }
        }
    }
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

    ned = S*(S-1)/2;
    nedm1 = ned - 1;

    init_subgraph_table();
    init_replicas();

    if (argc == 4) load_state(argv[3]);

    run();
    
    free_subgraph_table();

    return EXIT_SUCCESS;
}
