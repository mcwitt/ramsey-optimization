/* Genetic algorithm based on Goldberg's "Simple Genetic Algorithm" */

#include <stdio.h>
#include "dsfmt.h"
#include "sga.h"

#define SGA_RANDOM() dsfmt_genrand_close_open(&dsfmt)
#define SGA_RND(l, u) (u-l)*SGA_RANDOM() + l

#define MAX(x, y) ((x) > (y)) ? x : y
#define MIN(x, y) ((x) > (y)) ? y : x

dsfmt_t dsfmt;

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

/*
 * Select an individual from the population with probability proportional to
 * its fitness. First determine partitions by calling `set_parts'.
 */
static int select(double parts[])
{
    return bisect(parts, SGA_RANDOM(), 0, SGA_MAXPOP);
}

/* compute partitions for roulette-wheel selection */
static void set_parts(SGA_indiv_t pop[], double parts[])
{
    double sum = 0.;
    int i;

    parts[0] = 0.;
    for (i = 1; i < SGA_MAXPOP; i++) sum += (parts[i] = pop[i].fitness);
    for (i = 1; i < SGA_MAXPOP; i++) parts[i] = parts[i-1] + parts[i]/sum;
}

/* cross 2 parent strings at specified crossing site,
 * place in 2 child strings */
static void crossover(int parent1[], int parent2[],
                      int child1[],  int child2[],
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
static void mutate(int chrom[], double pmutate, int *nmutation)
{
    int i;

    for (i = 0; i < SGA_CHROMLEN; i++)
    {
        if (SGA_RANDOM() < pmutate)
        {
            chrom[i] = !chrom[i];
            *nmutation += 1;
        }
    }
}

static void update(SGA_indiv_t *p, int parent1, int parent2, int xsite)
{
    p->fitness = SGA_objfunc(p->chrom);
    p->parent1 = parent1;
    p->parent2 = parent2;
    p->xsite = xsite;
}

void SGA_init(uint32_t seed)
{
    /* init random number generator */
    dsfmt_init_gen_rand(&dsfmt, seed);
}

void SGA_init_pop(SGA_indiv_t pop[])
{
    SGA_indiv_t *p;
    int iind, j;

    for (iind = 0; iind < SGA_MAXPOP; iind++)
    {
        p = &pop[iind];
        
        /* initialize chromosome with random bits */
        for (j = 0; j < SGA_CHROMLEN; j++)
            p->chrom[j] = (SGA_RANDOM() < 0.5) ? 0 : 1;

        /* compute fitness; set parents and crossing site to default values */
        update(p, -1, -1, 0);
    }
}

void SGA_advance(SGA_indiv_t oldpop[], SGA_indiv_t newpop[],
                 SGA_params_t *params, SGA_stats_t *stats)
{
    double parts[SGA_MAXPOP];
    double f1, f2;
    int iind;
    int mate1, mate2, xsite = -1;

    stats->fitness_avg = 0.;
    stats->fitness_var = 0.;
    stats->fitness_min = 1e10;
    stats->fitness_max = -1e10;
    stats->ncross = 0;
    stats->nmutation = 0;

    set_parts(oldpop, parts);

    for (iind = 0; iind < params->npop; iind += 2)
    {
        mate1 = select(parts);
        mate2 = select(parts);

        if (SGA_RANDOM() < params->pcross)
        {
            xsite = SGA_RND(1, SGA_CHROMLEN-1);
            crossover(oldpop[mate1].chrom, oldpop[mate2 ].chrom,
                      newpop[iind ].chrom, newpop[iind+1].chrom,
                      xsite);
            stats->ncross += 1;
        }

        mutate(newpop[iind  ].chrom, params->pmutate, &stats->nmutation);
        mutate(newpop[iind+1].chrom, params->pmutate, &stats->nmutation);

        update(&newpop[iind  ], mate1, mate2, xsite);
        update(&newpop[iind+1], mate1, mate2, xsite);

        f1 = newpop[iind  ].fitness;
        f2 = newpop[iind+1].fitness;
        stats->fitness_avg += (f1 + f2);
        stats->fitness_var += (f1*f1 + f2*f2);
        stats->fitness_max = MAX(stats->fitness_min, MAX(f1, f2));
        stats->fitness_min = MIN(stats->fitness_min, MIN(f1, f2));
    }

    stats->fitness_avg /= params->npop;
    stats->fitness_var = stats->fitness_var / params->npop
                         - stats->fitness_avg * stats->fitness_avg;
}
