/* 
 * Energy annealing schedule for use with demon2.c
 * Energies decrease in a geometric sequence 
 */

#include "esched.h"

double R_er[NSTAGE_MAX][NEDR+1];
double R_es[NSTAGE_MAX][NEDS+1];

static void init_esched(double *e[NSTAGE_MAX], int nedt, int nstage);

void R_init_esched(int nstage)
{
    int i, n;
    double dx, x = 0.;

    er[0][0] = es[0][0] = 1.;
    for (i = 0; i < NEDR; i++) er[0][i] = es[0][i] = 0.;
    for (; i < NEDS; i++) es[0][i] = 0.;

    dx = 1./nstage;

    for (i = 0; i < nstage; i++)
    {
        x += dx;
        er[i][0] = es[i][0] = 1.;
        for (n = 1; n < NEDR; n++) e[n] = x*e[n-1];
        for (n = 1; n < NEDR; n++) e[n] = x*e[n-1];
    }
}

int main()
{
    R_init_esched(10);
    return 0;
}
