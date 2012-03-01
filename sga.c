/* Genetic algorithm based on Goldberg's "Simple Genetic Algorithm" */

#include "dSFMT.h"
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
static int select_mate(double parts[], int popsize)
{
    return bisect(parts, SGA_RANDOM(), 0, popsize);
}

/* compute partitions for roulette-wheel selection */
static void set_parts(SGA_indiv_t pop[], double parts[], int popsize)
{
    double sum = 0.;
    int i;

    parts[0] = 0.;
    for (i = 1; i < popsize; i++) sum += (parts[i] = pop[i].fitness);
    for (i = 1; i < popsize; i++) parts[i] = parts[i-1] + parts[i]/sum;
}

/* cross 2 parent strings at specified crossing site,
 * place in 2 child strings */
static void crossover(int parent1[], int parent2[],
                      int child1[],  int child2[],
                      int lchrom, int xsite)
{
    int i;

    for (i = 0; i < xsite; i++)
    {
        child1[i] = parent2[i];
        child2[i] = parent1[i];
    }

    for (; i < lchrom; i++)
    {
        child1[i] = parent1[i];
        child2[i] = parent2[i];
    }
}

/* flip each bit in a chromosome with probability pmutate */
static void mutate(int chrom[], int lchrom, double pmutate, int *nmutation)
{
    int i;

    for (i = 0; i < lchrom; i++)
    {
        if (SGA_RANDOM() < pmutate)
        {
            chrom[i] = !chrom[i];
            *nmutation += 1;
        }
    }
}

/* compute an individual's fitness and set parentage data */
static void init_indiv(SGA_indiv_t *p, int parent1, int parent2, int xsite)
{
    p->fitness = SGA_objfunc(p->chrom);
    p->parent1 = parent1;
    p->parent2 = parent2;
    p->xsite   = xsite;
}

void SGA_init(uint32_t seed)
{
    /* init random number generator */
    dsfmt_init_gen_rand(&dsfmt, seed);
}

void SGA_init_pop(SGA_indiv_t pop[], SGA_params_t *params)
{
    SGA_indiv_t *p;
    int iind, j;

    for (iind = 0; iind < params->popsize; iind++)
    {
        p = &pop[iind];
        
        /* initialize chromosome with random bits */
        for (j = 0; j < params->lchrom; j++)
            p->chrom[j] = (SGA_RANDOM() < 0.5) ? 0 : 1;

        /* compute fitness; set parents and crossing site to default values */
        init_indiv(p, -1, -1, 0);
    }
}

void SGA_advance(SGA_indiv_t oldpop[], SGA_indiv_t newpop[],
                 SGA_params_t *p, SGA_stats_t *s)
{
    double parts[SGA_MAXPOPSIZE];
    double f1, f2;
    int iind, mate1, mate2, xsite = 0;

    /* reset stats */
    s->fittest     =  NULL;
    s->fitness_avg =  0.;
    s->fitness_var =  0.;
    s->fitness_min =  1e10;
    s->fitness_max = -1e10;
    s->ncross      =  0;
    s->nmutation   =  0;

    set_parts(oldpop, parts, p->popsize);

    for (iind = 0; iind < p->popsize; iind += 2)
    {
        /* select mates with probability proportional to fitness */
        mate1 = select_mate(parts, p->popsize);
        mate2 = select_mate(parts, p->popsize);

        /* do crossover with probability pcross */
        if (SGA_RANDOM() < p->pcross)
        {
            xsite = SGA_RND(1, p->lchrom-1);
            s->ncross += 1;
        }
        else xsite = 0; /* copy without crossover */

        crossover(oldpop[mate1].chrom, oldpop[mate2 ].chrom,
                  newpop[iind ].chrom, newpop[iind+1].chrom,
                  p->lchrom, xsite);

        mutate(newpop[iind  ].chrom, p->lchrom, p->pmutate, &s->nmutation);
        mutate(newpop[iind+1].chrom, p->lchrom, p->pmutate, &s->nmutation);

        /* compute fitness and set parentage data for children */
        init_indiv(&newpop[iind  ], mate1, mate2, xsite);
        init_indiv(&newpop[iind+1], mate1, mate2, xsite);

        /* update stats */
        f1 = newpop[iind  ].fitness;
        f2 = newpop[iind+1].fitness;
        s->fitness_avg += (f1 + f2);
        s->fitness_var += (f1*f1 + f2*f2);
        s->fitness_min = MIN(s->fitness_min, MIN(f1, f2));
        if (f1 > s->fitness_max) { s->fitness_max = f1; s->fittest = &newpop[iind  ]; }
        if (f2 > s->fitness_max) { s->fitness_max = f2; s->fittest = &newpop[iind+1]; }
    }

    s->fitness_avg /= p->popsize;
    s->fitness_var = s->fitness_var / p->popsize
                         - s->fitness_avg * s->fitness_avg;
}
