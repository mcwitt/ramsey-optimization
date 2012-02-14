#include <stdio.h>
#include "ramsey2.h"

dsfmt_t R_rstate;

/* subs[ei] (subr[ei]) lists the NSGFES (NSGFER) complete S-subgraphs
 * (R-subgraphs) that include edge ei */
int *subr[NED];
int *subs[NED];

double R_er[NEDR];
double R_es[NEDS];
 
double der[NEDR-1]; /* der[i] = R_er[i] - R_er[i+1] */
double des[NEDS-1];

void init_tabs(int *sub[], int t, int nsgfe);
void free_tabs(int *sub[]);

void R_init(uint32_t seed)
{
    int i;

    /* init random number generator */
    dsfmt_init_gen_rand(&R_rstate, seed);

    for (i = 0; i < NEDR-1; i++) der[i] = R_er[i] - R_er[i+1];
    for (i = 0; i < NEDS-1; i++) des[i] = R_es[i] - R_es[i+1];

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
    for (j = 0; j < NSGR; j++) p->nb[j] = NEDR;
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
            p->en += R_flip_energy(p, j);
            p->sp[j] = -1;
            R_update(p, j);
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
            p->en += R_flip_energy(p, j);
            p->sp[j] *= -1;
            R_update(p, j);
        }
    }
}

double R_flip_energy(rep_t *p, int edge)
{
    double en = 0.;
    int isub;
    int *nr = p->nr;
    int *nb = p->nb;

    if (p->sp[edge] == 1)
    {
        for (isub = 0; isub < NSGFER; isub++)
            en += der[nb[subr[edge][isub]]-1];
        for (isub = 0; isub < NSGFES; isub++)
            en -= des[nr[subs[edge][isub]]];
    }
    else
    {
        for (isub = 0; isub < NSGFER; isub++)
            en -= der[nb[subr[edge][isub]]];
        for (isub = 0; isub < NSGFES; isub++)
            en += des[nr[subs[edge][isub]]-1];
    }

    return en;
}

void R_update(rep_t *p, int edge)
{
    int isub;
    int sp = p->sp[edge];

    for (isub = 0; isub < NSGFER; isub++) p->nb[subr[edge][isub]] += sp;
    for (isub = 0; isub < NSGFES; isub++) p->nr[subs[edge][isub]] -= sp;
}

void R_update_energy(rep_t *p)
{
    int j;

    for (j = 0; j < NEDR-1; j++) der[j] = R_er[j] - R_er[j+1];
    for (j = 0; j < NEDS-1; j++) des[j] = R_es[j] - R_es[j+1];

    p->en = (double) NSGS;
    for (j = 0; j < NSGR; j++) p->nb[j] = NEDR;
    for (j = 0; j < NSGS; j++) p->nr[j] = 0;

    for (j = 0; j < NED; j++) if (p->sp[j] == -1)
    {
        p->sp[j] = 1;   /* hack */
        p->en += R_flip_energy(p, j);
        p->sp[j] = -1;   /* hack */
        R_update(p, j);
    }
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

