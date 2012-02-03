/* Undefined constants (computed by compile script)
 *      NV      : number of vertices
 *      R       : red clique size
 *      S       : blue clique size
 *      NED     : number of edges (=NV(NV-1)/2)
 *      NSGR    : number of subgraphs with R vertices (=binomial(NV, S))
 *      NSGS    : number of subgraphs with S vertices (=binomial(NV, S))
 *      NSGFER  : number of subgraphs with R vertices including a given edge
 *      NSGFES  : number of subgraphs with S vertices including a given edge
 *                  (=binomial(NV-2, S-2))
 */

#include "dSFMT.h"

#if !(  defined(NV)     && defined(R)       && defined(S) &&\
        defined(NED)    && defined(NSGR)    && defined(NSGS) &&\
        defined(NSGFER) && defined(NSGFES))
#error Missing required definitions. See dSFMT.h.
#endif

#define R_RAND() dsfmt_genrand_close_open(&rng_state)
        
/* Strucure to store the configuration of one replica */
typedef struct
{
    int sp[NED];
    int h2[NED];
    int nbr[NSGR];  /* number of blue edges in each R-subgraph */
    int nbs[NSGS];  /* number of blue edges in each S-subgraph */
    int en;         /* number of blue S-cliques and red R-cliques */
} rep_t;

extern dsfmt_t rng_state;   /* state of random number generator (RNG) */

/* Call R_init() first, and R_finalize() to free memory when done */
void R_init(uint32_t seed);
void R_finalize();

/* Initialize replica with all edges blue */
void R_init_replica(rep_t *p);  

/*
 * Initialize replica in a random configuration with equal numbers of red and
 * blue edges on average
 */
void R_init_replica_random(rep_t *p);

/*
 * Initialize replica using configuration from a graph file. If the file
 * specifies a graph with fewer than NV vertices, initialize the unspecified
 * edges randomly with equal probabilities for red and blue. Returns the number
 * of edges that were initialized with random values.
 */
int R_init_replica_from_file(rep_t *p, char filename[]);

/* Randomizes spins with indices less than imask */
void R_randomize(rep_t *p, int imask);

/*
 * Update local fields after a spin flip (called AFTER spin assigned new value)
 *      p   : pointer to affected replica
 *      ei  : index of flipped spin
 */
void R_update_fields(rep_t *p, int ei);

/* Save graph to file */
void R_save_graph(int sp[], char filename[]);