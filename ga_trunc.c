/* Genetic algorithm based on Goldberg's "Simple Genetic Algorithm" but using
 * Truncation Selection (i.e. pick the k best individuals for reproduction)
 * instead of Fitness Proportional Probability Selection.
 */

#include <assert.h>
#include <math.h>
#include "dSFMT.h"
#include "ga_trunc.h"

#ifdef CROSS2
#warning Using 2-point crossover
#endif

#define RANDOM() dsfmt_genrand_close_open(&dsfmt)
#define RND(l, u) (u-l)*RANDOM() + l

dsfmt_t dsfmt;

/* comparison function used by rank */
int cmp_index(void *thunk, const void *a, const void *b) {
    double *arr = (double*) thunk;
    int ia = *(int*)a, ib = *(int*)b;
    return (arr[ia] - arr[ib] < 0.) ? 1 : -1;   /* largest first */
}

/* rank fitnesses and create a table of indices by rank */
static void rank(int index[], double fitness[], int n)
{
    int i;

    i = n; while (--i) index[i] = i; index[0] = 0;
    qsort_r(index, n, sizeof(int), fitness, cmp_index);
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
}

void GA_init(GA_t *ga, int popsize, int lchrom,
              double (*objfunc)(GA_allele_t*), double (*fitfunc)(double),
              int k, double pcross, double pmutate, uint32_t seed)
{
    int i, j;

    assert(popsize % 2 == 0);

    ga->popsize = popsize;
    ga->lchrom  = lchrom;
    ga->objfunc = objfunc;
    ga->fitfunc = fitfunc;
    ga->k = k;
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
    int index[GA_MAXPOPSIZE];
    int i, r, mate1, mate2;

    *ncross    = 0;
    *nmutation = 0;

    rank(index, ga->fitness, ga->popsize);

    /* create new generation using crossover and mutation */
    for (i = 0; i < ga->popsize; i += 2)
    {
        /* select mates from k fittest individuals */
        r = (int) RND(0, ga->k); mate1 = index[r];
        r = (int) RND(0, ga->k); mate2 = index[r];

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
