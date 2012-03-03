/* Genetic algorithm based on Goldberg's "Simple Genetic Algorithm" */

#include <math.h>
#include "dSFMT.h"
#include "sga.h"

#define FMULT 1.5
#define C     3.0

#define MAX(a, b) ((a) > (b)) ? a : b

#define SGA_RANDOM() dsfmt_genrand_close_open(&dsfmt)
#define SGA_RND(l, u) (u-l)*SGA_RANDOM() + l

dsfmt_t dsfmt;

/*
 * "sigma truncation" (Goldberg 124)--translate fitness values such that the
 * average is c*sigma, then set negative values to zero
 */
static void sigmatrunc(double in[], double out[], int len, double c,
                       double *avg, double *var, double *min, double *max)
{
    double shift, x;
    int i;

    shift = c*sqrt(*var) - *avg;

    *var = 0.;
    *avg += shift;
    *max += shift;
    if ((*min += shift) < 0.) *min = 0.;

    for (i = 0; i < len; i++)
    {
        x = in[i] + shift;
        if (x < 0.) { *avg = *avg - x/len; x = 0.; }
        else *var += x*x;
        out[i] = x;
    }

    *var = *var/len - (*avg) * (*avg);
}

/* do linear scaling of fitnesses as described in Goldberg pp. 78-79 */
static void scalepop(double in[], double out[], int len, double mult,
                     double *avg, double *var, double *min, double *max)
{
    double delta, a, b;
    int i;

    sigmatrunc(in, out, len, C, avg, var, min, max);

    /* determine linear scaling coefficients */
    /* non-negative test */
    if (*min > (FMULT * (*avg) - (*max)) / (FMULT - 1.))
    {
        /* normal scaling */
        delta = *max - *avg;
        a = (FMULT - 1.) * (*avg) / delta;
        b = *avg * ((*max) - FMULT*(*avg)) / delta;
    }
    else
    {
        /* scale as much as possible */
        delta = *avg - *min;
        a = *avg / delta;
        b = -(*avg) * (*min) / delta;
    }

    /* apply scaling and calculate new variance*/
    for (i = 0; i < len; i++) out[i] = a*out[i] + b;
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
    return bisect(parts, SGA_RANDOM(), 0, popsize);
}

/* compute partitions for roulette-wheel selection */
static void preselect(double fitness[], double parts[],
                      double sumfitness, int popsize)
{
    int i;

    parts[0] = fitness[0] / sumfitness;
    for (i = 1; i < popsize; i++) parts[i] = parts[i-1] + fitness[i]/sumfitness;
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

void SGA_init(SGA_t *sga, int popsize, int lchrom,
              double pcross, double pmutate, uint32_t seed)
{
    double f;
    int i, j;

    sga->popsize = popsize;
    sga->lchrom  = lchrom;
    sga->pcross  = pcross;
    sga->pmutate = pmutate;

    /* init random number generator */
    dsfmt_init_gen_rand(&dsfmt, seed);

    sga->chrom = sga->chrom1;
    sga->next = sga->chrom2;

    /* create random population */
    sga->fittest =  0;
    sga->favg    =  0.;
    sga->fvar    =  0.;
    sga->fmin    =  1e10;
    sga->fmax    = -1e10;

    for (i = 0; i < popsize; i++)
    {
        /* initialize chromosome with random bits */
        for (j = 0; j < lchrom; j++)
            sga->chrom[i][j] = (SGA_RANDOM() < 0.5) ? 0 : 1;

        f = sga->fitness[i] = sga->objective[i] = SGA_objfunc(sga->chrom[i]);
        sga->favg += f;
        sga->fvar += f*f;
        if (f > sga->fmax) {sga->fmax = f; sga->fittest = i; }
        if (f < sga->fmin) {sga->fmin = f; }
    }

    sga->favg /= popsize;
    sga->fvar = sga->fvar / popsize - (sga->favg)*(sga->favg);

    scalepop(sga->fitness, sga->fitness, popsize, FMULT,
             &sga->favg, &sga->fvar, &sga->fmin, &sga->fmax);

}

void SGA_advance(SGA_t *sga, int *ncross, int *nmutation)
{
    double f;
    double parts[SGA_MAXPOPSIZE];
    int i, mate1, mate2, xsite = 0;
    int (*swap)[SGA_MAXPOPSIZE];

    *ncross    =  0;
    *nmutation =  0;

    preselect(sga->fitness, parts, sga->popsize * sga->favg, sga->popsize);

    /* create new generation using crossover and mutation */
    for (i = 0; i < sga->popsize; i += 2)
    {
        /* select mates with probability proportional to fitness */
        mate1 = select_mate(parts, sga->popsize);
        mate2 = select_mate(parts, sga->popsize);

        /* do crossover with probability pcross */
        if (SGA_RANDOM() < sga->pcross)
        {
            xsite = SGA_RND(1, sga->lchrom-1);
            *ncross += 1;
        }
        else xsite = sga->lchrom; /* copy without crossover */

        crossover(sga->chrom[mate1], sga->chrom[mate2],
                  sga->next[i], sga->next[i+1],
                  sga->lchrom, xsite);

        mutate(sga->next[i  ], sga->lchrom, sga->pmutate, nmutation);
        mutate(sga->next[i+1], sga->lchrom, sga->pmutate, nmutation);
    }

    /* swap pointers so that sga->chrom points to the new generation */
    swap = sga->chrom;
    sga->chrom = sga->next;
    sga->next = swap;

    /* compute fitness and statistics for new generation */
    sga->favg  =  0.;
    sga->fvar  =  0.;
    sga->fmax  = -1e10;
    sga->fmin  =  1e10;

    for (i = 0; i < sga->popsize; i++)
    {
        f = sga->fitness[i] = sga->objective[i] = SGA_objfunc(sga->chrom[i]);
        sga->favg += f;
        sga->fvar += f*f;
        if (f > sga->fmax) {sga->fmax = f; sga->fittest = i; }
        if (f < sga->fmin) {sga->fmin = f; }
    }

    sga->favg /= sga->popsize;
    sga->fvar = sga->fvar / sga->popsize - (sga->favg)*(sga->favg);

    scalepop(sga->fitness, sga->fitness, sga->popsize, FMULT,
             &sga->favg, &sga->fvar, &sga->fmin, &sga->fmax);
}
