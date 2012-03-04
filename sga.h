#include "stdint.h"

#define SGA_MAXPOPSIZE  1000
#define SGA_MAXLCHROM   1000

typedef int allele_t;

typedef struct
{
    /* buffers to store chromosome data */
    allele_t b1[SGA_MAXPOPSIZE][SGA_MAXLCHROM];
    allele_t b2[SGA_MAXPOPSIZE][SGA_MAXLCHROM];

    /* pointers to chromosome data for current and next generation */
    allele_t (*chrom)[SGA_MAXPOPSIZE], (*nextg)[SGA_MAXPOPSIZE];

    double (*objfunc)(allele_t*);

    double objective[SGA_MAXPOPSIZE];   /* raw value of objective function */
    double fitness[SGA_MAXPOPSIZE];     /* scaled fitness */

    int popsize, lchrom;
    double pcross, pmutate;
    int fittest;
    double fmin, fmax, favg, fvar;
} SGA_t;

/* create a population of random individuals */
void SGA_init(SGA_t *sga, int popsize, int lchrom,
              double (*objfunc)(allele_t*),
              double pcross, double pmutate, uint32_t seed);

/* advance one generation */
void SGA_advance(SGA_t *sga, int *ncross, int *nmutation);
