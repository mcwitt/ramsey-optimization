#include "stdint.h"

#define SGA_MAXPOPSIZE  1000
#define SGA_MAXLCHROM   1000

/* application-specific function to be defined in external file */
double SGA_objfunc(int chrom[]);

typedef struct
{
    int chrom1[SGA_MAXPOPSIZE][SGA_MAXLCHROM];
    int chrom2[SGA_MAXPOPSIZE][SGA_MAXLCHROM];
    int (*chrom)[SGA_MAXPOPSIZE], (*next)[SGA_MAXPOPSIZE];
    double objective[SGA_MAXPOPSIZE];
    double fitness[SGA_MAXPOPSIZE];
    int popsize, lchrom;
    double pcross, pmutate;
    int fittest;
    double fmin, fmax, favg, fvar;
} SGA_t;

/* create a population of random individuals */
void SGA_init(SGA_t *sga, int popsize, int lchrom,
              double pcross, double pmutate, uint32_t seed);

/* advance one generation */
void SGA_advance(SGA_t *sga, int *ncross, int *nmutation);
