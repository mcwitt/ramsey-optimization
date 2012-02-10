#include <stdio.h>
#include "ramsey.h"

dsfmt_t rng_state;

/* subs[ei] (subr[ei]) lists the NSGFES (NSGFER) complete S-subgraphs
 * (R-subgraphs) that include edge ei */
int *subr[NED];
int *subs[NED];
 
/* edgs[si] (edgr[si]) lists the edges of S-subgraph (R-subgraph) si */
int *edgr[NSGR];
int *edgs[NSGS];

int nedr = R*(R-1)/2;       /* number of edges in an R-subgraph */
int neds = S*(S-1)/2;       /* number of edges in an S-subgraph */

void init_tabs(int *sub[], int *edg[], int t, int nsgfe, int nedrs);
void free_tabs(int *sub[], int *edg[], int nsg);

void R_init(uint32_t seed)
{
    /* init random number generator */
    dsfmt_init_gen_rand(&rng_state, seed);

    init_tabs(subr, edgr, R, NSGFER, nedr);
    init_tabs(subs, edgs, S, NSGFES, neds);
}

void R_finalize()
{
    free_tabs(subr, edgr, NSGR);
    free_tabs(subs, edgs, NSGS);
}

void R_init_replica(rep_t *p)
{
    int j;

    p->en = NSGS;
    for (j = 0; j < NSGR; j++) p->nb[j] = nedr;
    for (j = 0; j < NSGS; j++) p->nr[j] = 0;

    for (j = 0; j < NED; j++)
    {
        p->sp[j] = 1;
        p->h2[j] = -NSGFES;
    }
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
            p->sp[j] = -1;
            p->en += p->h2[j];
            R_update_fields(p, j);
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
            p->en += p->sp[j]*p->h2[j];
            p->sp[j] *= -1;
            R_update_fields(p, j);
        }
    }
}

void R_update_fields(rep_t *p, int ei)
{
    int si, j, ej, n;
    int *sp = p->sp;
    int *h2 = p->h2;
    int *nr = p->nr;
    int *nb = p->nb;

    if (sp[ei] == 1)
    {
        for (si = 0; si < NSGFER; si++)
        {
            n = nb[subr[ei][si]] += 1;

            if (n > 2) continue;

            if (n == 2) /* destroyed an incomplete red clique */
            {
                for (j = 0; j < nedr; j++)
                {
                    ej = edgr[subr[ei][si]][j];
                    if (sp[ej] == 1 && ej != ei) { h2[ej] -= 1; break; }
                }
            }
            else        /* (n = 1) destroyed a red clique */
            {
                h2[ei] += 1;
                for (j = 0; j < nedr; j++) h2[edgr[subr[ei][si]][j]] -= 1;
            }
        }

        for (si = 0; si < NSGFES; si++)
        {
            n = nr[subs[ei][si]] -= 1;

            if (n > 1) continue;

            if (n == 1) /* created an incomplete blue clique */
            {
                for (j = 0; j < neds; j++)
                {
                    ej = edgs[subs[ei][si]][j];
                    if (sp[ej] == -1) { h2[ej] -= 1; break; }
                }
            }
            else        /* (n = 0) created a blue clique */
            {
                h2[ei] += 1;
                for (j = 0; j < neds; j++) h2[edgs[subs[ei][si]][j]] -= 1;
            }
        }
    }
    else
    {
        for (si = 0; si < NSGFER; si++)
        {
            n = nb[subr[ei][si]] -= 1;

            if (n > 1) continue;

            if (n == 1) /* created an incomplete red clique */
            {
                for (j = 0; j < nedr; j++)
                {
                    ej = edgr[subr[ei][si]][j];
                    if (sp[ej] == 1) { h2[ej] += 1; break; }
                }
            }
            else        /* (n = 0) created a red clique */
            {
                h2[ei] -= 1;
                for (j = 0; j < nedr; j++) h2[edgr[subr[ei][si]][j]] += 1;
            }
        }

        for (si = 0; si < NSGFES; si++)
        {
            n = nr[subs[ei][si]] += 1;

            if (n > 2) continue;

            if (n == 2) /* destroyed an incomplete blue clique */
            {
                for (j = 0; j < neds; j++)
                {
                    ej = edgs[subs[ei][si]][j];
                    if (sp[ej] == -1 && ej != ei) { h2[ej] += 1; break; }
                }
            }
            else        /* (n = 1) destroyed a blue clique */
            {
                h2[ei] -= 1;
                for (j = 0; j < neds; j++) h2[edgs[subs[ei][si]][j]] += 1;
            }
        }
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

void init_tabs(int *sub[], int *edg[], int t, int nsgfe, int nedrs)
{
    int ps[NED];    /* current positions in subgraph arrays */
    int pe;         /* current position in edge array */
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

        edg[si] = (int*) malloc(nedrs * sizeof(int));
        pe = 0;

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
                edg[si][pe++] = ei;
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

void free_tabs(int *sub[], int *edg[], int nsg)
{
    int i;

    for (i = 0; i < NED; i++)
        free(sub[i]);
    for (i = 0; i < nsg; i++)
        free(edg[i]);
}

