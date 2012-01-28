/* 
 * File: ramsey.h
 *
 * Version 4
 *
 * Author: Matt Wittmann <mwittman@ucsc.edu>
 *
 * Description: Parallel tempering Monte Carlo code which attempts to minimize
 * the number of r-cliques and s-independent sets of a graph given that it has
 * N_v vertices. Energy is defined to be the sum of the number of r-cliques and
 * s-independent sets. If for a given input (r, s, N_v) we find a zero-energy
 * state, this implies that R(r, s) > N_v.
 *
 * Type definitions and global variable and function declarations
 */

#include <time.h>
#include "dSFMT.h"

#define MAX_NT          32
#define MAX_SWEEPS      90000
#define WRITE_INTERVAL  90000

#define URAND() dsfmt_genrand_close_open(&rstate)

/* The following data structure stores the state of a single replica */
typedef struct
{
    int sp[NED];    /* "spin" of each edge [+1 (blue) or -1 (red)] */
    int h2[NED];    /* (doubled) local field for each spin */
    int nbr[NSGR];  /* number of blue edges in each R-subgraph */
    int nbs[NSGS];  /* number of blue edges in each S-subgraph */
    int en;         /* number of blue S-cliques and red R-cliques */
} rep_t;

extern int nedr;    /* number of edges in an R-subgraph */
extern int neds;    /* number of edges in an S-subgraph */
extern int nedsm1;  /* neds minus one (stored for efficiency) */
extern int nedsm2;  /* neds minus two */

/*
 * subs[ei] (subr[ei]) lists the complete S-subgraphs
 * (R-subgraphs) that include edge ei
 */
extern int *subr[NED];
extern int *subs[NED];

/* edgs[si] (edgr[si]) lists the edges of S-subgraph (R-subgraph) si */
extern int *edgr[NSGR];
extern int *edgs[NSGS];

extern rep_t reps[MAX_NT]; /* array of parallel tempering (PT) replicas */
extern int ri[MAX_NT];     /* replica indices in order of increasing temperature */

extern int nsweeps; /* number of sweeps */
extern int min;     /* lowest energy found */

extern int nt;                 /* number of PT copies */
extern int nswaps[MAX_NT];     /* number of swaps between each pair of temperatures */
extern double T[MAX_NT];       /* array of temperatures */
extern double mbeta[MAX_NT];   /* negative inverse temperatures */

extern dsfmt_t rstate; /* state of random number generator (RNG) */
extern uint32_t rseed; /* seed used to initialize RNG */

#ifndef NOTIME
extern clock_t start;  /* start time */
#endif

void init_reps();   /* initialize each replica in a random state */

/* set up look-up tables `sub' and `edg' */
void init_tabs(int *sub[], int *edg[], int t, int nsg, int nsgfe, int nedsg);

void run();     /* start the simulation loop */

void sweep();   /* update each spin once */

void temper();  /* attempt PT swaps */

/* call to update local fields after spin `ei' is flipped */
void update_fields(int ei, int sp[], int nbr[], int nbs[], int h2[]);

/* input/output functions defined in `io.c' */
void save_graph(int sp[NED], char filename[]);
void save_state(char filename[]);
void load_state(char filename[]);
void print_status();

