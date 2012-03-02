/* Interface for energy annealing schedules for use with demon2.c */

#include "ramsey2.c"

#define NSTAGE_MAX 100

extern double R_er[NSTAGE_MAX][NEDR+1];
extern double R_es[NSTAGE_MAX][NEDS+1];

void R_init_esched(int nstage);

