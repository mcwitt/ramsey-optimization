/* Genetic algorithm based on Goldberg's "Simple Genetic Algorithm" */

#include "dSFMT.h"
#include "sga.h"

#define SGA_RANDOM() dsfmt_genrand_close_open(&dsfmt)
#define SGA_RND(l, u) (u-l)*SGA_RANDOM() + l

#define MAX(x, y) ((x) > (y)) ? x : y
#define MIN(x, y) ((x) > (y)) ? y : x

dsfmt_t dsfmt;

/* do linear scaling of fitnesses as described in Goldberg pp. 78-79 */
static void linscale(SGA_indiv_t pop[], SGA_stats_t *st, SGA_params_t *pm)
{
    double delta, a, b;
    double umin, umax, uavg, fmul, f;

    umin = st->minfitness;
    umax = st->maxfitness;
    uavg = st->sumfitness / pm->popsize;
    fmul = SGA_FMULT;

    /* DETERMINE LINEAR SCALING COEFFICIENTS */
    /* non-negative test */
    if (umin > (fmul*uavg - umax) / (fmul - 1.))
    {
        /* normal scaling */
        delta = umax - uavg;
        a = (fmul - 1.) * uavg / delta;
        b = uavg * (umax - fmul*uavg) / delta;
    }
    else
    {
        /* scale as much as possible */
        delta = uavg - umin;
        a = uavg / delta;
        b = -uavg * umin / delta;
    }

    /* APPLY SCALING */
    st->maxfitness = a * st->maxfitness + b;
    st->minfitness = a * st->minfitness + b;
    st->sumfitness  = 0.;
    st->sumfitness2 = 0.;
    for (i = 0; i < pm->popsize; i++)
    {
        f = a*pop[i].fitness + b;
        pop[i].fitness = f;
        st->sumfitness += f;
        st->sumfitness2 += f*f;
    }
}

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
 * its fitness. First determine partitions by calling `preselect'.
 */
static int select_mate(double parts[], int popsize)
{
    return bisect(parts, SGA_RANDOM(), 0, popsize);
}

/* compute partitions for roulette-wheel selection */
static void preselect(SGA_indiv_t pop[], double parts[],
                      double sumfitness, int popsize)
{
    int i;

    parts[0] = pop[0].fitness / sumfitness;
    for (i = 1; i < popsize; i++) parts[i] = parts[i-1] + pop[i].fitness/sumfitness;
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
        child1[i] = parent1[i];
        child2[i] = parent2[i];
    }

    for (; i < lchrom; i++)
    {
        child1[i] = parent2[i];
        child2[i] = parent1[i];
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
    p->objective = SGA_objfunc(p->chrom);
    p->parent1 = parent1;
    p->parent2 = parent2;
    p->xsite   = xsite;
}

void SGA_init(uint32_t seed)
{
    /* init random number generator */
    dsfmt_init_gen_rand(&dsfmt, seed);
}

void SGA_init_pop(SGA_indiv_t pop[], SGA_stats_t *st, SGA_params_t *pm)
{
    double f;
    int i, j;

    st->fittest     =  0;
    st->sumfitness  =  0.;
    st->sumfitness2 =  0.;
    st->maxfitness  =  1e10;
    st->minfitness  = -1e10;

    for (i = 0; i < pm->popsize; i++)
    {
        /* initialize chromosome with random bits */
        for (j = 0; j < pm->lchrom; j++)
            pop[i].chrom[j] = (SGA_RANDOM() < 0.5) ? 0 : 1;

        /* compute fitness; set parents and crossing site to default values */
        init_indiv(&pop[i], -1, -1, 0);

        f = pop[i].fitness;
        st->sumfitness  += f;
        st->sumfitness2 += f*f;
        st->minfitness = MIN(st->minfitness, f);
        if (f > st->maxfitness)
            { st->fittest = i; st->maxfitness = f; }
    }
}

void SGA_advance(SGA_indiv_t oldpop[], SGA_indiv_t newpop[],
                 SGA_stats_t *st, SGA_params_t *pm)
{
    double f1, f2;
    double parts[SGA_MAXPOPSIZE];
    int i, mate1, mate2, xsite = 0;

    preselect(oldpop, parts, st->sumfitness, pm->popsize);

    st->fittest     =  0;
    st->ncross      =  0;
    st->nmutation   =  0;
    st->sumfitness  =  0.;
    st->sumfitness2 =  0.;
    st->maxfitness  = -1e10;
    st->minfitness  =  1e10;

    for (i = 0; i < pm->popsize; i += 2)
    {
        /* select mates with probability proportional to fitness */
        mate1 = select_mate(parts, pm->popsize);
        mate2 = select_mate(parts, pm->popsize);

        /* do crossover with probability pcross */
        if (SGA_RANDOM() < pm->pcross)
        {
            xsite = SGA_RND(1, pm->lchrom-1);
            st->ncross += 1;
        }
        else xsite = pm->lchrom; /* copy without crossover */

        crossover(oldpop[mate1].chrom, oldpop[mate2 ].chrom,
                  newpop[i    ].chrom, newpop[i+1   ].chrom,
                  pm->lchrom, xsite);

        mutate(newpop[i  ].chrom, pm->lchrom, pm->pmutate, &st->nmutation);
        mutate(newpop[i+1].chrom, pm->lchrom, pm->pmutate, &st->nmutation);

        /* compute fitness and set parentage data for children */
        init_indiv(&newpop[i  ], mate1, mate2, xsite);
        init_indiv(&newpop[i+1], mate1, mate2, xsite);

        /* accumulate stats */
        f1 = newpop[i  ].fitness;
        f2 = newpop[i+1].fitness;
        st->sumfitness  += (f1 + f2);
        st->sumfitness2 += (f1*f1 + f2*f2);
        st->minfitness = MIN(st->minfitness, MIN(f1, f2));

        /* check for new fittest individual */
        if (f1 > st->maxfitness)
            { st->maxfitness = f1; st->fittest = i;   }
        if (f2 > st->maxfitness)
            { st->maxfitness = f2; st->fittest = i+1; }
    }
}
