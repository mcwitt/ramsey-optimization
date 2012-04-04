#include "dSFMT.h"
#include "defs.h"   /* constant definitions generated by gendefs.py */

#if (R > S)
#error Must have R < S
#endif

/* Strucure to store the configuration of one replica */
typedef struct
{
    int sp[NED];
    int nb[NSGR];       /* number of blue edges in each R-subgraph */
    int nr[NSGS];       /* number of red edges in each S-subgraph */
    double h2[NED];     /* local field */
    double der[NEDR];   /* der[i] = er[i] - er[i+1] */
    double des[NEDS];
    double en;          /* number of blue S-cliques and red R-cliques */
} R_replica_t;

/* Call R_init() first, and R_finalize() to free memory when done */
void R_init(uint32_t seed);
void R_finalize();

/* Initialize replica with all edges blue */
void R_init_replica(R_replica_t *p);  

/*
 * Initialize replica using configuration from a .graph file. If the file
 * specifies a graph with fewer than NV vertices, the unspecified edges will be
 * blue. Returns the number of spins read from the file (this number can be
 * passed as the "mask" argument to R_randomize to prevent these spins from
 * being randomized).
 */
int R_init_replica_from_file(R_replica_t *p, char filename[]);

/* Randomize spins with indices greater than "mask" */
void R_randomize(R_replica_t *p, double p_red, int mask);

/* Flip spin with index iedge (note this does NOT update the energy) */
void R_flip(R_replica_t *p, int iedge);

/* Compute the energy of a graph from scratch */
double R_energy(int sp[NED], double er[NEDR+1], double es[NEDS+1]);

/* 
 * Set subgraph energy levels. er[nb] the energy of an R-subgraph with nb blue
 * edges. Note that this changes the total energy in general.
 */
void R_set_energies(R_replica_t *p, double er[NEDR+1], double es[NEDS+1]);

/* Save graph to file */
void R_save_graph(int sp[], char filename[]);
