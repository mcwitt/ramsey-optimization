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

#if (R > S)
#error Please change definitions so that R < S.
#endif

#define NEDR R*(R-1)/2  /* number of edges in an R-subgraph */
#define NEDS S*(S-1)/2  /* number of edges in an S-subgraph */
#define R_RAND() dsfmt_genrand_close_open(&R_rstate)

/* Strucure to store the configuration of one replica */
typedef struct
{
    int sp[NED];
    int nb[NSGR];       /* number of blue edges in each R-subgraph */
    int nr[NSGS];       /* number of red edges in each S-subgraph */
    double der[NEDR];   /* der[i] = er[i] - er[i+1] */
    double des[NEDS];
    double en;          /* number of blue S-cliques and red R-cliques */
} rep_t;

extern dsfmt_t R_rstate;    /* state of random number generator (RNG) */

/* Call R_init() first, and R_finalize() to free memory when done */
void R_init(uint32_t seed);
void R_finalize();

/* Initialize replica with all edges blue */
void R_init_replica(rep_t *p);  

/*
 * Initialize replica using configuration from a graph file. If the file
 * specifies a graph with fewer than NV vertices, the unspecified edges will be
 * blue. Returns the number of unspecified edges (this number can be passed as
 * imask to R_randomize to randomize just the unspecified edges).
 */
int R_init_replica_from_file(rep_t *p, char filename[]);

/* Randomize spins with indices less than imask */
void R_randomize(rep_t *p, double p_red, int imask);

/* Compute the energy to flip a spin */
double R_flip_energy(rep_t *p, int edge);

/* Must be called after each single spin flip to update state variables */
void R_update(rep_t *p, int edge);

/* 
 * Set subgraph energy levels. er is an array of energy levels for the
 * R-subgraphs, er[i] being the energy with i blue edges.
 */
void R_set_energies(rep_t *p, double er[NEDR+1], double es[NEDS+1]);

/* Save graph to file */
void R_save_graph(int sp[], char filename[]);
