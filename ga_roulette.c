/* Genetic algorithm based on Goldberg's "Simple Genetic Algorithm" */

#include <assert.h>
#include <math.h>
#include "dSFMT.h"
#include "ga_roulette.h"

#define FMULT   2.0     /* linear scaling parameter */
#define C       2.0     /* sigma truncation parameter */
#define DELTA   0.0     /* sigma truncation parameter */
#define EPSILON 10e-15  /* set to roughly machine precision */

#ifdef NOSCALE
#warning Linear scaling of fitnesses disabled
#endif

#ifdef NOTRUNC
#warning Sigma truncation disabled
#endif

#ifdef CROSS2
#warning Using 2-point crossover
#endif

#define RANDOM() dsfmt_genrand_close_open(&dsfmt)
#define RND(l, u) (u-l)*RANDOM() + l

dsfmt_t dsfmt;

#ifndef NOTRUNC
/*
 * "Sigma truncation" (Goldberg 124)--translate fitness values to make the
 * average c*sigma + delta, then set negative values to zero.
 * Nonzero delta avoids issues when fitness values are all equal
 */
static void sigmatrunc(double x[], int len, double c, double delta,
                       double *avg, double *var, double *min, double *max)
{
    double shift, y;
    int i;

    shift = c*sqrt(*var) - *avg + delta;

    *avg = *var = 0.;
    *max += shift;
    if ((*min += shift) < 0.) *min = 0.;

    for (i = 0; i < len; i++)
    {
        if ((y = x[i] + shift) < 0.) y = 0.;
        else { *avg += y; *var += y*y; }
        x[i] = y;
    }

    *avg /= len;
    *var = *var/len - (*avg)*(*avg);
}
#endif

#ifndef NOSCALE
/* do linear scaling of fitnesses as described in Goldberg pp. 78-79 */
static void linscale(double x[], int len, double mult,
                     double *avg, double *var, double *min, double *max)
{
    double delta, a, b;
    int i;

    /* determine linear scaling coefficients */
    /* non-negative test */
    if (*min > (mult * (*avg) - (*max)) / (mult - 1.))
    {
        /* normal scaling */
        delta = *max - *avg;
        a = (mult - 1.) * (*avg) / delta;
        b = *avg * ((*max) - mult*(*avg)) / delta;
    }
    else
    {
        /* scale as much as possible */
        delta = *avg - *min;
        a = *avg / delta;
        b = -(*avg) * (*min) / delta;
    }

    /* apply scaling and calculate new variance (average is unchanged) */
    for (i = 0; i < len; i++) x[i] = a*x[i] + b;
    *max = a * (*max) + b;
    *min = a * (*min) + b;
    *var = (a*a - 1)*(*avg)*(*avg) + 2*a*b*(*avg) + a*a*(*var) + b*b;
}
#endif

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
    return bisect(parts, RANDOM(), 0, popsize);
}

/* compute partitions for roulette-wheel selection */
static void preselect(double fitness[], double parts[],
                      double sumfitness, int popsize)
{
    int i;

    parts[0] = fitness[0] / sumfitness;
    for (i = 1; i < popsize; i++) parts[i] = parts[i-1] + fitness[i]/sumfitness;
}

/*
 * do 2-point crossover of 2 parent strings using specified crossing sites
 * x1 < x2, place results in 2 child strings
 */
static void cross(GA_allele_t parent1[], GA_allele_t parent2[],
                  GA_allele_t child1[],  GA_allele_t child2[],
                  int lchrom, int x1, int x2)
{
    int i;

    for (i = 0; i < x1; i++)
    {
        child1[i] = parent1[i];
        child2[i] = parent2[i];
    }
    
    for (; i < x2; i++)
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

/* 
 * with probability pcross, do 1-point crossover of 2 parent strings using
 * random crossing site, place results in 2 child strings; otherwise copy
 * parent strings
 */
#ifndef CROSS2
static void crossover(GA_allele_t parent1[], GA_allele_t parent2[],
                      GA_allele_t child1[],  GA_allele_t child2[],
                      int lchrom, double pcross, int *ncross)
{
    int xsite;

    if (RANDOM() < pcross) { xsite = (int) RND(1, lchrom-1); *ncross += 1; }
    else xsite = lchrom;    /* copy without crossover */
    cross(parent1, parent2, child1, child2, lchrom, xsite, lchrom);
}

/* 
 * with probability pcross, do 2-point crossover of 2 parent strings using
 * random crossing sites, place results in 2 child strings; otherwise copy
 * parent strings
 */
#else
static void crossover2(GA_allele_t parent1[], GA_allele_t parent2[],
                       GA_allele_t child1[],  GA_allele_t child2[],
                       int lchrom, double pcross, int *ncross)
{
    int x1, x2, swap;

    if (RANDOM() < pcross)
    {
        x1 = (int) RND(1, lchrom-1);
        x2 = (int) RND(1, lchrom-1);
        if (x1 > x2) { swap = x1; x1 = x2; x2 = swap; }
        *ncross += 1;
    }
    else x1 = x2 = lchrom;
    cross(parent1, parent2, child1, child2, lchrom, x1, x2);
}
#endif

/* flip each bit in a chromosome with probability pmutate */
static void mutate(GA_allele_t chrom[], int lchrom,
                   double pmutate, int *nmutation)
{
    int i;

    for (i = 0; i < lchrom; i++)
    {
        if (RANDOM() < pmutate)
        {
            chrom[i] = !chrom[i];
            *nmutation += 1;
        }
    }
}

/* calculate fitnesses and statistics */
void update_stats(double (*objfunc)(GA_allele_t*), double (*fitfunc)(double),
                  GA_allele_t (*chrom)[GA_MAXPOPSIZE],
                  double objective[], double fitness[], int popsize,
                  double *avg, double *var, double *min, double *max,
                  int *fittest)
{
    double f;
    int i;

    *avg =  0.;
    *var =  0.;
    *min =  1e10;
    *max = -1e10;

    for (i = 0; i < popsize; i++)
    {
        objective[i] = objfunc(chrom[i]);
        f = fitness[i] = fitfunc(objective[i]);
        *avg += f;
        *var += f*f;
        if (f > *max) { *max = f; *fittest = i; }
        if (f < *min) *min = f;
    }

    *avg /= popsize;
    *var = *var / popsize - (*avg)*(*avg);

#ifndef NOTRUNC
    sigmatrunc(fitness, popsize, C, DELTA, avg, var, min, max);
#endif
#ifndef NOSCALE
    if (*var > EPSILON) linscale(fitness, popsize, FMULT, avg, var, min, max);
#endif
}

void GA_init(GA_t *ga, int popsize, int lchrom,
              double (*objfunc)(GA_allele_t*), double (*fitfunc)(double),
              double pcross, double pmutate, uint32_t seed)
{
    int i, j;

    assert(popsize % 2 == 0);

    ga->popsize = popsize;
    ga->lchrom  = lchrom;
    ga->objfunc = objfunc;
    ga->fitfunc = fitfunc;
    ga->pcross  = pcross;
    ga->pmutate = pmutate;

    /* init random number generator */
    dsfmt_init_gen_rand(&dsfmt, seed);

    ga->chrom = ga->b1;
    ga->nextg = ga->b2;

    /* create random population */
    for (i = 0; i < popsize; i++)
        for (j = 0; j < lchrom; j++)
            ga->chrom[i][j] = (RANDOM() < 0.5) ? 0 : 1;

    update_stats(
        ga->objfunc, ga->fitfunc,
        ga->chrom, ga->objective, ga->fitness, ga->popsize,
        &ga->favg, &ga->fvar, &ga->fmin, &ga->fmax, &ga->fittest
    );
}

void GA_advance(GA_t *ga, int *ncross, int *nmutation)
{
    GA_allele_t (*swap)[GA_MAXPOPSIZE];
    double parts[GA_MAXPOPSIZE];
    int i, mate1, mate2;

    *ncross    = 0;
    *nmutation = 0;

    preselect(ga->fitness, parts, ga->popsize * ga->favg, ga->popsize);

    /* create new generation using crossover and mutation */
    for (i = 0; i < ga->popsize; i += 2)
    {
        /* select mates with probability proportional to fitness */
        mate1 = select_mate(parts, ga->popsize);
        mate2 = select_mate(parts, ga->popsize);

        /* do crossover with probability pcross */
#ifndef CROSS2
        crossover(ga->chrom[mate1], ga->chrom[mate2],
                  ga->nextg[i], ga->nextg[i+1],
                  ga->lchrom, ga->pcross, ncross);
#else
        crossover2(ga->chrom[mate1], ga->chrom[mate2],
                   ga->nextg[i], ga->nextg[i+1],
                   ga->lchrom, ga->pcross, ncross);
#endif

        mutate(ga->nextg[i  ], ga->lchrom, ga->pmutate, nmutation);
        mutate(ga->nextg[i+1], ga->lchrom, ga->pmutate, nmutation);
    }

    /* swap pointers so that ga->chrom points to the new generation */
    swap = ga->chrom;
    ga->chrom = ga->nextg;
    ga->nextg = swap;

    /* compute fitness and statistics for new generation */
    update_stats(
        ga->objfunc, ga->fitfunc,
        ga->chrom, ga->objective, ga->fitness, ga->popsize,
        &ga->favg, &ga->fvar, &ga->fmin, &ga->fmax, &ga->fittest
    );
}
