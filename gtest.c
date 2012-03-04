#include <stdio.h>
#include <math.h>
#include "sga.h"

#define LCHROM  30
#define POPSIZE 30

double SGA_objfunc(int chrom[])
{
    int i;
    double powerof2 = 1., x = 0., coef = 1073741823.;

    for (i = 0; i < LCHROM; i++)
    {
        if (chrom[i] == 1) x += powerof2;
        powerof2 *= 2.;
    }

    return pow(x/coef, 10);
}

void chrom2str(int chrom[], char str[])
{
    int i;

    for (i = 0; i < LCHROM; i++)
        str[i] = (chrom[i] == 1) ? '1' : '0';
}

int main()
{
    SGA_t sga;
    int ncross, nmutation;

    SGA_init(&sga, POPSIZE, 30, 0.6, 0.00333, 123);
    SGA_advance(&sga, &ncross, &nmutation);
}
