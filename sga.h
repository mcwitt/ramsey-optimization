#include "stdint.h"

#define SGA_MAXPOPSIZE  1000
#define SGA_MAXLCHROM   1000

/* application-specific function to be defined in external file */
double SGA_objfunc(int chrom[]);

typedef struct
{
    int chrom[SGA_MAXLCHROM];
    double fitness;
    int parent1, parent2, xsite;
} SGA_indiv_t;

typedef struct
{
    int popsize;
    int lchrom;
    double pcross;
    double pmutate;
} SGA_params_t;

typedef struct
{
    SGA_indiv_t *fittest;   /* pointer to fittest individual */
    double fitness_avg, fitness_var, fitness_min, fitness_max;
    int ncross;
    int nmutation;
} SGA_stats_t;

/* call first to initialize */
void SGA_init(uint32_t seed);

/* create a population of random individuals */
void SGA_init_pop(SGA_indiv_t pop[], SGA_params_t *params);

/* advance oldpop by one generation; store result in newpop */
void SGA_advance(SGA_indiv_t oldpop[], SGA_indiv_t newpop[],
                 SGA_params_t *params, SGA_stats_t *stats);