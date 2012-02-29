#define SGA_POPSIZE 100

typedef struct
{
    allele_t chrom[SGA_CHROMLEN];
    pheno_t pheno;
    double fitness;
    int parent1, parent2, xsite;
} SGA_indiv_t;

typedef struct
{
    double pcross;
    double pmutate;
} SGA_params_t;

typedef struct
{
    int ncross;
    int nmutation;
} SGA_stats_t;

/* create a population of random individuals */
void SGA_init_pop(SGA_indiv_t pop[]);

/* advance oldpop by one generation; store result in newpop */
void SGA_advance(SGA_indiv_t oldpop[], SGA_indiv_t newpop[],
                 SGA_params_t *params, SGA_stats_t *stats);
