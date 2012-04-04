#include <stdio.h>
#include "dSFMT.h"
#include "ramsey2h.h"

#define RANDOM() dsfmt_genrand_close_open(&dsfmt)
dsfmt_t dsfmt;

/* subs[iedge] (subr[iedge]) lists the NSGFES (NSGFER) complete S-subgraphs
 * (R-subgraphs) that include edge iedge */
int *subr[NED];
int *subs[NED];
 
/* edgs[isub] (edgr[isub]) lists the edges of S-subgraph (R-subgraph) isub */
int *edgr[NSGR];
int *edgs[NSGS];

int nedr = R*(R-1)/2;       /* number of edges in an R-subgraph */
int neds = S*(S-1)/2;       /* number of edges in an S-subgraph */

static void init_tabs(int *sub[], int *edg[], int t, int nsgfe, int nedt)
{
    int nsub[NED];  /* number of subgraphs processed for each edge */
    int nedg;       /* number of edges of the current subgraph processed */
    int v[S+2];     /* vertices of the current subgraph */
    int iedg, isub; /* edge and subgraph indices */
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
     * Art of Computer Programming Vol. 4, Fasc. 3 (all-caps labels correspond
     * to labels in the book)
     */

    /* INITIALIZE */
    isub = 0;
    v[t] = NV;
    v[t+1] = 0;
    for (j = 0; j < t; j++) v[j] = j;

    for (;;)
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
void R_init(uint32_t seed)
{
    /* init random number generator */
    dsfmt_init_gen_rand(&dsfmt, seed);

    init_tabs(subr, edgr, R, NSGFER, nedr);
    init_tabs(subs, edgs, S, NSGFES, neds);
}

void R_finalize()
{
    int i;

    for (i = 0; i < NED; i++)  { free(subr[i]); free(subs[i]); }
    for (i = 0; i < NSGR; i++) { free(edgr[i]); free(edgs[i]); }
    for (; i < NSGS; i++) free(edgs[i]);
}

void R_init_replica(R_replica_t *p)
{
    int j;

    for (j = 0; j < NSGR; j++) p->nb[j] = nedr;
    for (j = 0; j < NSGS; j++) p->nr[j] = 0;

    for (j = 0; j < NED; j++)
    {
        p->sp[j] = 1;
        p->h2[j] = (double) -NSGFES;
    }

    p->der[0] = p->des[0] = 1.;
    for (j = 1; j < NEDR; j++) p->der[j] = 0.;
    for (j = 1; j < NEDS; j++) p->des[j] = 0.;

    p->en = (double) NSGS;
}

int R_init_replica_from_file(R_replica_t *p, char filename[])
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

void R_randomize(R_replica_t *p, double p_red, int mask)
{
    int j;

    for (j = mask; j < NED; j++)
    {
        if (RANDOM() < p_red)
        {
            p->en += p->sp[j]*p->h2[j];
            R_flip(p, j);
        }
    }
}

void R_flip(R_replica_t *p, int edge)
{
    double dhr, dhb;
    double *h2 = p->h2;
    int isub, sg, iedge;
    int *edges;
    int *sp = p->sp;
    int *nr = p->nr;
    int *nb = p->nb;

    if ((sp[edge] *= -1) == 1)
    {
        for (isub = 0; isub < NSGFER; isub++)
        {
            sg = subr[edge][isub];
            nb[sg] += 1;
            dhr = -p->der[nb[sg]-1];
            /* dhr = R_er[nb[sg]] - R_er[nb[sg]-1]; */
            dhb = -p->der[nb[sg]-2];
            /* dhb = R_er[nb[sg]-1] - R_er[nb[sg]-2]; */
            edges = edgr[subr[edge][isub]];

            for (iedge = 0; iedge < nedr; iedge++) if (edges[iedge] != edge)
                h2[edges[iedge]] += ((sp[edges[iedge]] == 1) ? dhb : dhr);
        }

        for (isub = 0; isub < NSGFES; isub++)
        {
            sg = subs[edge][isub];
            nr[sg] -= 1;
            dhr = des[nr[sg]-1];
            /* dhr = R_es[nr[sg]-1] - R_es[nr[sg]]; */
            dhb = des[nr[sg]];
            /* dhb = R_es[nr[sg]] - R_es[nr[sg]+1]; */
            edges = edgs[subs[edge][isub]];

            for (iedge = 0; iedge < neds; iedge++) if (edges[iedge] != edge)
                h2[edges[iedge]] -= ((sp[edges[iedge]] == 1) ? dhb : dhr);
        }
    }
    else
    {
        for (isub = 0; isub < NSGFER; isub++)
        {
            sg = subr[edge][isub];
            nb[sg] -= 1;
            dhr = der[nb[sg]];
            /* dhr = R_er[nb[sg]] - R_er[nb[sg]+1]; */
            dhb = der[nb[sg]-1];
            /* dhb = R_er[nb[sg]-1] - R_er[nb[sg]]; */
            edges = edgr[subr[edge][isub]];

            for (iedge = 0; iedge < nedr; iedge++) if (edges[iedge] != edge)
                h2[edges[iedge]] += ((sp[edges[iedge]] == 1) ? dhb : dhr);
        }

        for (isub = 0; isub < NSGFES; isub++)
        {
            sg = subs[edge][isub];
            nr[sg] += 1;
            dhr = -des[nr[sg]-2];
            /* dhr = R_es[nr[sg]-1] - R_es[nr[sg]-2]; */
            dhb = -des[nr[sg]-1];
            /* dhb = R_es[nr[sg]] - R_es[nr[sg]-1]; */
            edges = edgs[subs[edge][isub]];

            for (iedge = 0; iedge < neds; iedge++) if (edges[iedge] != edge)
                h2[edges[iedge]] -= ((sp[edges[iedge]] == 1) ? dhb : dhr);
        }
    }
}

void R_set_energies(R_replica_t *p, double er[], double es[])
{
    int j;

    for (j = 0; j < NEDR; j++) p->der[j] = er[j] - er[j+1];
    for (j = 0; j < NEDS; j++) p->des[j] = es[j] - es[j+1];

    p->en = (double) NSGS;
    for (j = 0; j < NSGR; j++) p->nb[j] = NEDR;
    for (j = 0; j < NSGS; j++) p->nr[j] = 0;

    for (j = 0; j < NED; j++) if (p->sp[j] == -1)
    {
        /* want energy relative to graph with all edges blue, so need to
         * temporarily set edge to blue to compute the correct energy */
        p->sp[j] = 1;
        p->en += R_flip_energy(p, j);
        p->sp[j] = -1;

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

double R_energy(int sp[NED], double er[NEDR+1], double es[NEDS+1])
{
    double energy = 0.;
    int isub, iedg, sum;

    for (isub = 0; isub < NSGR; isub++)
    {
        sum = 0;
        for (iedg = 0; iedg < NEDR; iedg++) sum += sp[edgr[isub][iedg]];
        energy += er[(NEDR+sum)>>1];
    }

    for (isub = 0; isub < NSGS; isub++)
    {
        sum = 0;
        for (iedg = 0; iedg < NEDS; iedg++) sum += sp[edgs[isub][iedg]];
        energy += es[(NEDS-sum)>>1];
    }

    return energy;
}
