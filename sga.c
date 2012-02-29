/* Genetic algorithm based on Goldberg's "Simple Genetic Algorithm" */

#include <stdio.h>
#include "sga.h"

/*
 * Select an individual from the population with probability proportional to
 * its fitness. First determine partitions by calling `set_parts'.
 */
#define SELECT(parts) bisect(parts, RANDOM(), 0, SGA_POPSIZE)

/* return the insertion point for x to maintain sorted order of a */
static int bisect(double a[], double x, int l, int r)
{
    int mid;

    while (r-l > 1)
    {
        mid = (l+r)/2;
        if (a[mid] < x) l = mid;
        else r = mid;
    }

    return (x < a[l]) ? l : r;
}

static int select(double parts[])
{
    return bisect(parts, RANDOM(), 0, SGA_POPSIZE);
}

/* compute partitions for roulette-wheel selection */
static void set_parts(SGA_indiv_t pop[], double parts[])
{
    double sum = 0.;
    int i;

    parts[0] = 0.;
    for (i = 1; i < SGA_POPSIZE; i++) sum += (parts[i] = pop[i].fitness);
    for (i = 1; i < SGA_POPSIZE; i++) parts[i] = parts[i-1] + parts[i]/sum;
}

/* cross 2 parent strings at specified crossing site,
 * place in 2 child strings */
static void crossover(allele_t parent1[], allele_t parent2[],
                      allele_t child1[],  allele_t child2[],
                      int xsite)
{
    int i;

    for (i = 0; i < xsite; i++)
    {
        child1[i] = parent1[i];
        child2[i] = parent2[i];
    }

    for (; i < SGA_CHROMLEN; i++)
    {
        child1[i] = parent2[i];
        child2[i] = parent1[i];
    }
}

/* apply mutations to chromosome with bit-flipping probability pmutate */
static void mutate(allele_t chrom[], double pmutate, int *nmutation)
{
    int i;

    for (i = 0; i < SGA_CHROMLEN; i++)
    {
        if (RANDOM() < pmutation)
        {
            chrom[i] = !chrom[i];
            *nmutation += 1;
        }
    }
}

static void update(SGA_indiv_t *p, int parent1, int parent2, int xsite)
{
    decode(p->chrom, &p->pheno);
    p->fitness = objfunc(&p->pheno);
    p->parent1 = parent1;
    p->parent2 = parent2;
    p->xsite = xsite;
}

void SGA_init_pop(SGA_indiv_t pop[])
{
    int iind, j;

    for (iind = 0; iind < SGA_POPSIZE; iind++)
    {
        p = &pop[iind];
        
        /* initialize chromosome with random bits */
        for (j = 0; j < SGA_CHROMLEN; j++)
            p->chrom[j] = (RANDOM() < 0.5) ? 0 : 1;

        /* compute phenotype and fitness */
        update(p, -1, -1, 0);
    }
}

void SGA_advance(SGA_indiv_t oldpop[], SGA_indiv_t newpop[],
                 SGA_params_t *params, SGA_stats_t *stats);
{
    double parts[SGA_POPSIZE];
    int iind;

    set_parts(oldpop, parts);

    for (iind = 0; iind < SGA_POPSIZE; iind += 2)
    {
        mate1 = select(parts);
        mate2 = select(parts);

        if (RANDOM() < params->pcross)
        {
            xsite = RND(1, SGA_CHROMLEN-1);
            crossover(oldpop[mate1].chrom, oldpop[mate2].chrom,
                      newpop[i    ].chrom, newpop[i + 1].chrom,
                      xsite);
            stats->ncross += 1;
        }

        mutate(newpop[i  ].chrom, params->pmutate, &stats->nmutation);
        mutate(newpop[i+1].chrom, params->pmutate, &stats->nmutation);

        update(&newpop[i  ], mate1, mate2, xsite);
        update(&newpop[i+1], mate1, mate2, xsite);
    }
}


int main(int argc, char *argv[])
{
    SGA_indiv_t oldpop[SGA_POPSIZE], newpop[SGA_POPSIZE], *pswap;
    SGA_params_t params;
    SGA_stats_t stats;
    double parts[];
    int mate1, mate2, xsite;
    int nmutation;
    int igen, iind;

    params->pcross = 0.6;
    params->pmutate = 0.001;

    stats->ncross = 0;
    stats->nmutation = 0;

    if (argc != 7 && argc != 8)
    {
        fprintf(stderr, "Usage: %s npop ngen pcross pmutate seed", argv[0]);
        fprintf(stderr, "Compiled for (%d, %d, %d)\n", R, S, NV);
        exit(EXIT_FAILURE);
    }

    init_pop(pop);

    for (igen = 0; igen < ngen; igen++)
    {
        SGA_advance(oldpop, newpop, &params, &stats);

        /* swap pointers to population arrays so that oldpop becomes newpop
         * (former oldpop will be overwritten in next iteration) */
        pswap = oldpop;
        oldpop = newpop;
        newpop = pswap;
    }

    return EXIT_SUCCESS;
}
