/*
 * Undefined constants
 *      NV  : number of vertices
 *      S   : clique size
 */

#define NT_MAX 24
#define WRITE_INTERVAL 10
#define OUTFILE "zero.data"

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
int ri[NT_MAX];     /* replica indices in order of ascending temperature */

int nt;                 /* number of PT copies */
int nswaps[NT_MAX];     /* number of swaps between each pair of temperatures */
double T[NT_MAX];       /* array of temperatures */
double mbeta[NT_MAX];   /* negative inverse temperatures */

dsfmt_t rstate; /* state of random number generator (RNG) */
uint32_t rseed; /* seed used to initialize RNG */

/* include debugging functions */
#ifdef DEBUG
#include "debug.c"
#endif

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

int flip_energy(int sp, int *sub, int *nb)
{
    int i, nbi, delta;

    delta = 0;

    if (sp == 1)
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

void update_nb(int sp, int *sub, int *nb)
{
    int i;

    for (i = 0; i < nsgfe; i++) nb[sub[i]] += sp;
}

/* initialize each replica with a random configuration */
void init_reps()
{
    rep_t *p;
    int it, j, k;

    for (it = 0; it < nt; it++)
    {
        ri[it] = it;
        p = &reps[it];
        p->energy = nsg;
        p->nb = (int*) malloc(nsg * sizeof(int));
        nswaps[it] = 0;

        for (j = 0; j < nsg; j++) p->nb[j] = ned;

        for (k = 0; k < NV; k++)
        {
            for (j = 0; j < k; j++)
            {
                if (URAND() > 0.5) p->sp[j][k] = 1;
                else
                {
                    p->sp[j][k] = -1;
                    p->energy += flip_energy(1, sub[j][k], p->nb);
                    update_nb(-1, sub[j][k], p->nb);
                }
            }
        }
    }
}

void free_reps()
{
    int it;

    for (it = 0; it < nt; it++)
        free(reps[it].nb);
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
                delta = flip_energy(p->sp[j][k], sub[j][k], p->nb);

                /* flip with Metropolis probability */
                if (delta <= 0 || URAND() < exp(mbeta[it]*delta))
                {
                    p->energy += delta;
                    update_nb(p->sp[j][k] *= -1, sub[j][k], p->nb);
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
    fprintf(fp, "%d\n", S);
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

void run()
{
    char filename[256];
    int it, nsweeps, done;

    nsweeps = 0;
    done = 0;

    while (! done && nsweeps < 50)
    {
        sweep();
        nsweeps++;

        temper();

        for (it=0; it<nt; it++)
            printf("%3d ", reps[ri[it]].energy);
        for (it=1; it<nt; it++)
            printf("%3.2f ", (float) nswaps[it]/nsweeps);
        printf("\n");

        if (nsweeps % WRITE_INTERVAL == 0)
        {
            sprintf(filename, "%d-%d-%d_%d_%d.bin",
                    NV, S, S, rseed, nsweeps/WRITE_INTERVAL);
            save_state(filename);
        }

        for (it = 0; it < nt; it++)
        {
            if (reps[ri[it]].energy == 0)
            {
                printf("Found zero-energy ground state!\n");
                printf("Graph saved to %s\n", OUTFILE);
                printf("N_sweeps = %d\n", nsweeps);
                save_graph(reps[ri[it]].sp, OUTFILE);
                done = 1;
                break;
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
    while (fscanf(infile, "%lf", &t) != EOF && nt < NT_MAX)
    {
        T[nt] = t;
        mbeta[nt] = -1./t;
        nt++;
    }
    assert(nt > 1);

    init_globals();
    init_reps();
    if (argc == 4) load_state(argv[3]);

    run();
    
    free_reps();
    free_globals();

    return EXIT_SUCCESS;
}
