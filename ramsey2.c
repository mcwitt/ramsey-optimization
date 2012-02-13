#include <stdio.h>
#include "ramsey2.h"

dsfmt_t rng_state;

/* subs[ei] (subr[ei]) lists the NSGFES (NSGFER) complete S-subgraphs
 * (R-subgraphs) that include edge ei */
int *subr[NED];
int *subs[NED];
 
int nedr = R*(R-1)/2;       /* number of edges in an R-subgraph */
int neds = S*(S-1)/2;       /* number of edges in an S-subgraph */

void init_tabs(int *sub[], int t, int nsgfe);
void free_tabs(int *sub[]);

void R_init(uint32_t seed)
{
    /* init random number generator */
    dsfmt_init_gen_rand(&rng_state, seed);

    init_tabs(subr, R, NSGFER);
    init_tabs(subs, S, NSGFES);
}

void R_finalize()
{
    free_tabs(subr);
    free_tabs(subs);
}

void R_init_replica(rep_t *p)
{
    int j;

    p->en = (double) NSGS;
    for (j = 0; j < NSGR; j++) p->nb[j] = nedr;
    for (j = 0; j < NSGS; j++) p->nr[j] = 0;

    for (j = 0; j < NED; j++) p->sp[j] = 1;
}

void R_init_replica_random(rep_t *p)
{
    R_init_replica(p);
    R_randomize(p, NED);
}

int R_init_replica_from_file(rep_t *p, char filename[])
{
    FILE *fp;
    int ned, sp, j, imask;

    R_init_replica(p);

    if (! (fp = fopen(filename, "r")))
    {
        fprintf(stderr, "Error opening file: %s\n", filename);
        exit(EXIT_FAILURE);
    }

    if (! fscanf(fp, "%d", &ned))
    {
        fprintf(stderr, "Error while reading file: %s\n", filename);
        exit(EXIT_FAILURE);
    }

    ned = ned*(ned-1)/2;
    imask = NED - ned;      /* number of spins unspecified by input */
    assert(imask >= 0);
    R_randomize(p, imask);  /* randomize unspecified spins */

    /* read remaining spins from input */
    j = imask;
    while (fscanf(fp, "%d", &sp) != EOF && j < NED)
    {
        if (sp == 0)
        {
            p->en += R_h2(p, j);
            p->sp[j] = -1;
        }

        j++;
    }

    fclose(fp);

    return imask;
}

void R_randomize(rep_t *p, int imask)
{
    int j;

    for (j = 0; j < imask; j++)
    {
        if (R_RAND() < 0.5)
        {
            p->en += p->sp[j]*R_h2(p, j);
            p->sp[j] *= -1;
        }
    }
}

double er[] = {1., 0, 0, 0, 0, 0};
double es[] = {1., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

double R_h2(rep_t *p, int edge)
{
    double h2 = 0.;
    int isub;
    int sub;
    int *nr = p->nr;
    int *nb = p->nb;

    for (isub = 0; isub < NSGFER; isub++)
    {
        sub = subr[edge][isub];
        h2 += er[nb[sub]-1] - er[nb[sub]];  /* energy to create red edge */
    }

    for (isub = 0; isub < NSGFES; isub++)
    {
        sub = subs[edge][isub];
        h2 += es[nr[sub]+1] - es[nr[sub]];  /* energy to destroy blue edge */
    }

    return h2;
}

void R_save_graph(int sp[NED], char filename[])
{
    FILE *fp;
    int i;

    fp = fopen(filename, "w");
    fprintf(fp, "%d\n", NV);
    
    for (i = 0; i < NED; i++)
        fprintf(fp, "%d\n", (sp[i] == 1) ? 1 : 0);

    fclose(fp);
}

void init_tabs(int *sub[], int t, int nsgfe)
{
    int ps[NED];    /* current positions in subgraph arrays */
    int c[t+2];     /* array of vertices of the current subgraph */
    int ei, si;     /* edge index, subgraph index */
    int j, k;

    for (j = 0; j < NED; j++)
    {
        sub[j] = (int*) malloc(nsgfe * sizeof(int));
        ps[j] = 0;
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
                 * add subgraph si to list for edge ei
                 * add edge ei to list for subgraph si
                 */

                sub[ei][ps[ei]++] = si;
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

void free_tabs(int *sub[])
{
    int i;

    for (i = 0; i < NED; i++)
        free(sub[i]);
}

