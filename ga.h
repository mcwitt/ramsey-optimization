/* Genetic algorithm based on Goldberg's "Simple Genetic Algorithm"
 *
 * Author: Matt Wittmann <mwittman@ucsc.edu>
 *
 * Optional compilation flags:
 *  NOSCALE : don't perform linear scaling of fitness values before selection
 *  NOTRUNC : don't use sigma truncation (will generate errors if objective
 *              function returns negative values
 *  CROSS2  : use 2-point crossover
 */

#include "stdint.h"

#define GA_MAXPOPSIZE  1000
#define GA_MAXLCHROM   1000

typedef char GA_allele_t;

typedef struct
{
    /* buffers to store chromosome data */
    GA_allele_t b1[GA_MAXPOPSIZE][GA_MAXLCHROM];
    GA_allele_t b2[GA_MAXPOPSIZE][GA_MAXLCHROM];

    /* pointers to chromosome data for current and next generation */
    GA_allele_t (*chrom)[GA_MAXPOPSIZE], (*nextg)[GA_MAXPOPSIZE];

    double (*objfunc)(GA_allele_t*);   /* pointer to objective function */
    double (*fitfunc)(double);          /* pointer to fitness function */

    double objective[GA_MAXPOPSIZE];   /* objective function values */
    double fitness[GA_MAXPOPSIZE];     /* scaled fitnesses */

    int popsize, lchrom;
    double pcross, pmutate;
    double fmin, fmax, favg, fvar;
    int fittest;
} GA_t;

/* create a population of random individuals */
void GA_init(GA_t *sga, int popsize, int lchrom,
              double (*objfunc)(GA_allele_t*), double (*fitfunc)(double),
              double pcross, double pmutate, uint32_t seed);

/* advance one generation */
void GA_advance(GA_t *sga, int *ncross, int *nmutation);
