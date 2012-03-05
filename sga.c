/* Genetic algorithm based on Goldberg's "Simple Genetic Algorithm" */

#include <assert.h>
#include <math.h>
#include "dSFMT.h"
#include "sga.h"

#define FMULT   2.0
#define C       2.0
#define DELTA   0.0
#define EPSILON 10e-15  /* set to roughly machine precision */
//#define SGA_C2

#define RANDOM() dsfmt_genrand_close_open(&dsfmt)
#define RND(l, u) (u-l)*RANDOM() + l

dsfmt_t dsfmt;

/*
 * "sigma truncation" (Goldberg 124)--translate fitness values to make the
 * average c*sigma + delta, then set negative values to zero
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

    /* apply scaling and calculate new variance*/
    for (i = 0; i < len; i++) x[i] = a*x[i] + b;
    *max = a * (*max) + b;
    *min = a * (*min) + b;
    *var = (a*a - 1)*(*avg)*(*avg) + 2*a*b*(*avg) + a*a*(*var) + b*b;
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
static void cross(SGA_allele_t parent1[], SGA_allele_t parent2[],
                  SGA_allele_t child1[],  SGA_allele_t child2[],
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
#ifndef SGA_C2
static void crossover(SGA_allele_t parent1[], SGA_allele_t parent2[],
                      SGA_allele_t child1[],  SGA_allele_t child2[],
                      int lchrom, double pcross, int *ncross)
{
    int xsite;

    if (RANDOM() < pcross) { xsite = RND(1, lchrom-1); *ncross += 1; }
    else xsite = lchrom;    /* copy without crossover */
    cross(parent1, parent2, child1, child2, lchrom, xsite, lchrom);
}

/* 
 * with probability pcross, do 2-point crossover of 2 parent strings using
 * random crossing sites, place results in 2 child strings; otherwise copy
 * parent strings
 */
#else
static void crossover2(SGA_allele_t parent1[], SGA_allele_t parent2[],
                       SGA_allele_t child1[],  SGA_allele_t child2[],
                       int lchrom, double pcross, int *ncross)
{
    int x1, x2, swap;

    if (RANDOM() < pcross)
    {
        x1 = RND(1, lchrom-1);
        x2 = RND(1, lchrom-1);
        if (x1 > x2) { swap = x1; x1 = x2; x2 = swap; }
        *ncross += 1;
    }
    else x1 = x2 = lchrom;
    cross(parent1, parent2, child1, child2, lchrom, x1, x2);
}
#endif

/* flip each bit in a chromosome with probability pmutate */
static void mutate(SGA_allele_t chrom[], int lchrom,
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
void update_stats(double (*objfunc)(SGA_allele_t*), double (*fitfunc)(double),
                  SGA_allele_t (*chrom)[SGA_MAXPOPSIZE],
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

    sigmatrunc(fitness, popsize, C, DELTA, avg, var, min, max);
    if (*var > EPSILON) linscale(fitness, popsize, FMULT, avg, var, min, max);
}

void SGA_init(SGA_t *sga, int popsize, int lchrom,
              double (*objfunc)(SGA_allele_t*), double (*fitfunc)(double),
              double pcross, double pmutate, uint32_t seed)
{
    int i, j;

    assert(popsize % 2 == 0);

    sga->popsize = popsize;
    sga->lchrom  = lchrom;
    sga->objfunc = objfunc;
    sga->fitfunc = fitfunc;
    sga->pcross  = pcross;
    sga->pmutate = pmutate;

    /* init random number generator */
    dsfmt_init_gen_rand(&dsfmt, seed);

    sga->chrom = sga->b1;
    sga->nextg = sga->b2;

    /* create random population */
    for (i = 0; i < popsize; i++)
        for (j = 0; j < lchrom; j++)
            sga->chrom[i][j] = (RANDOM() < 0.5) ? 0 : 1;

    update_stats(
        sga->objfunc, sga->fitfunc,
        sga->chrom, sga->objective, sga->fitness, sga->popsize,
        &sga->favg, &sga->fvar, &sga->fmin, &sga->fmax, &sga->fittest
    );
}

void SGA_advance(SGA_t *sga, int *ncross, int *nmutation)
{
    SGA_allele_t (*swap)[SGA_MAXPOPSIZE];
    double parts[SGA_MAXPOPSIZE];
    int i, mate1, mate2;

    *ncross    = 0;
    *nmutation = 0;

    preselect(sga->fitness, parts, sga->popsize * sga->favg, sga->popsize);

    /* create new generation using crossover and mutation */
    for (i = 0; i < sga->popsize; i += 2)
    {
        /* select mates with probability proportional to fitness */
        mate1 = select_mate(parts, sga->popsize);
        mate2 = select_mate(parts, sga->popsize);

        /* do crossover with probability pcross */
#ifndef SGA_C2
        crossover(sga->chrom[mate1], sga->chrom[mate2],
                  sga->nextg[i], sga->nextg[i+1],
                  sga->lchrom, sga->pcross, ncross);
#else
        crossover2(sga->chrom[mate1], sga->chrom[mate2],
                   sga->nextg[i], sga->nextg[i+1],
                   sga->lchrom, sga->pcross, ncross);
#endif

        mutate(sga->nextg[i  ], sga->lchrom, sga->pmutate, nmutation);
        mutate(sga->nextg[i+1], sga->lchrom, sga->pmutate, nmutation);
    }

    /* swap pointers so that sga->chrom points to the new generation */
    swap = sga->chrom;
    sga->chrom = sga->nextg;
    sga->nextg = swap;

    /* compute fitness and statistics for new generation */
    update_stats(
        sga->objfunc, sga->fitfunc,
        sga->chrom, sga->objective, sga->fitness, sga->popsize,
        &sga->favg, &sga->fvar, &sga->fmin, &sga->fmax, &sga->fittest
    );
}
