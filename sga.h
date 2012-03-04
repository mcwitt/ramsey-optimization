#include "stdint.h"

#define SGA_MAXPOPSIZE  1000
#define SGA_MAXLCHROM   1000

typedef char SGA_allele_t;

typedef struct
{
    /* buffers to store chromosome data */
    SGA_allele_t b1[SGA_MAXPOPSIZE][SGA_MAXLCHROM];
    SGA_allele_t b2[SGA_MAXPOPSIZE][SGA_MAXLCHROM];

    /* pointers to chromosome data for current and next generation */
    SGA_allele_t (*chrom)[SGA_MAXPOPSIZE], (*nextg)[SGA_MAXPOPSIZE];

    double (*objfunc)(SGA_allele_t*);   /* pointer to objective function */
    double (*fitfunc)(double);          /* pointer to fitness function */

    double objective[SGA_MAXPOPSIZE];   /* objective function values */
    double fitness[SGA_MAXPOPSIZE];     /* scaled fitnesses */

    int popsize, lchrom;
    double pcross, pmutate;
    double fmin, fmax, favg, fvar;
    int fittest;
} SGA_t;

/* create a population of random individuals */
void SGA_init(SGA_t *sga, int popsize, int lchrom,
              double (*objfunc)(SGA_allele_t*), double (*fitfunc)(double),
              double pcross, double pmutate, uint32_t seed);

/* advance one generation */
void SGA_advance(SGA_t *sga, int *ncross, int *nmutation);
