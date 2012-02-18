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

int R_init_replica_from_file(rep_t *p, char filename[])
{
    FILE *fp;
    int sp, j, ned;

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

    j = 0;

    while ((fscanf(fp, "%d", &sp) != EOF) && (j < NED))
    {
        if (sp == 0)
        {
            p->en += p->h2[j];
            R_flip(p, j);
        }

        j++;
    }

    fclose(fp);
    return j;
}

void R_randomize(rep_t *p, double p_red, int mask)
{
    int j;

    for (j = mask; j < NED; j++)
    {
        if (R_RAND() < p_red)
        {
            p->en += p->sp[j]*p->h2[j];
            R_flip(p, j);
        }
    }
}

void R_flip(rep_t *p, int ei)
{
    int si, j, n;
    int *sp = p->sp;
    int *h2 = p->h2;
    int *nr = p->nr;
    int *nb = p->nb;
    int *e;

    if ((sp[ei] *= -1) == 1)
    {
        for (si = 0; si < NSGFER; si++)
        {
            n = nb[subr[ei][si]] += 1;
            if (n > 2) continue;
            e = edgr[subr[ei][si]];

            if (n == 2) /* destroyed an incomplete red clique */
            {
                for (j = 0; j < nedr; j++)
                    if (sp[e[j]] == 1 && e[j] != ei) { h2[e[j]] -= 1; break; }
            }
            else        /* (n = 1) destroyed a red clique */
            {
                h2[ei] += 1;
                for (j = 0; j < nedr; j++) h2[e[j]] -= 1;
            }
        }

        for (si = 0; si < NSGFES; si++)
        {
            n = nr[subs[ei][si]] -= 1;
            if (n > 1) continue;
            e = edgs[subs[ei][si]];

            if (n == 1) /* created an incomplete blue clique */
            {
                for (j = 0; j < neds; j++)
                    if (sp[e[j]] == -1) { h2[e[j]] -= 1; break; }
            }
            else        /* (n = 0) created a blue clique */
            {
                h2[ei] += 1;
                for (j = 0; j < neds; j++) h2[e[j]] -= 1;
            }
        }
    }
    else
    {
        for (si = 0; si < NSGFER; si++)
        {
            n = nb[subr[ei][si]] -= 1;
            if (n > 1) continue;
            e = edgr[subr[ei][si]];

            if (n == 1) /* created an incomplete red clique */
            {
                for (j = 0; j < nedr; j++)
                    if (sp[e[j]] == 1) { h2[e[j]] += 1; break; }
            }
            else        /* (n = 0) created a red clique */
            {
                h2[ei] -= 1;
                for (j = 0; j < nedr; j++) h2[e[j]] += 1;
            }
        }

        for (si = 0; si < NSGFES; si++)
        {
            n = nr[subs[ei][si]] += 1;
            if (n > 2) continue;
            e = edgs[subs[ei][si]];

            if (n == 2) /* destroyed an incomplete blue clique */
            {
                for (j = 0; j < neds; j++)
                    if (sp[e[j]] == -1 && e[j] != ei) { h2[e[j]] += 1; break; }
            }
            else        /* (n = 1) destroyed a blue clique */
            {
                h2[ei] -= 1;
                for (j = 0; j < neds; j++) h2[e[j]] += 1;
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

void init_tabs(int *sub[], int *edg[], int t, int nsgfe, int nedt)
{
    int nsub[NED];      /* number of subgraphs processed for each edge */
    int nedg;          /* number of edges of the current subgraph processed */
    int v[S+2];         /* vertices of the current subgraph */
    int iedg, isub;    /* edge index, subgraph index */
    int j, k;

    for (j = 0; j < NED; j++)
    {
        sub[j] = (int*) malloc(nsgfe * sizeof(int));
        nsub[j] = 0;
    }

    /* 
     * Iterate over all subgraphs with t vertices (i.e. combinations of t
     * vertices)
     *
     * Algorithm to generate combinations adapted from Algorithm L in Knuth's
     * Art of Computer Programming Vol. 4, Fasc. 3 (all-caps labels
     * correspond to labels in the book)
     */

    /* INITIALIZE */
    isub = 0;
    v[t] = NV;
    v[t+1] = 0;
    for (j = 0; j < t; j++) v[j] = j;

    while (1)
    {
        /*
         * VISIT subgraph v_1 v_2 ... v_t
         * (algorithm guarantees that v_1 < v_2 < ... < v_t)
         */

        edg[isub] = (int*) malloc(nedt * sizeof(int));
        nedg = 0;

        /* iterate over edges in this subgraph */
        for (j = 0; j < t; j++)
        {
            for (k = 0; k < j; k++)
            {
                iedg = v[j]*(v[j]-1)/2 + v[k];
                sub[iedg][nsub[iedg]++] = isub; /* append isub to sub[iedg] */
                edg[isub][nedg++] = iedg;       /* append iedg to edg[isub] */
            }
        }

        isub++;   /* increment subgraph label */

        /* FIND j */
        j = 0;
        while (v[j] + 1 == v[j+1]) { v[j] = j; j++; }

        /* DONE? */
        if (j == t) break;

        v[j]++;
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

