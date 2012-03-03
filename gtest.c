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

void report(FILE *fp, SGA_indiv_t pop[])
{
    SGA_indiv_t *p;
    int i;
    char cstr[LCHROM];

    fprintf(fp, "%4s %8s %6s %*s %8s\n",
            "#", "parents", "xsite", LCHROM+2, "string", "fitness");

    for (i = 0; i < POPSIZE; i++)
    {
        p = &pop[i];
        chrom2str(p->chrom, cstr);
        fprintf(fp, "%4d (%2d, %2d) %6d %*s %8.4f\n",
                i,
                p->parent1,
                p->parent2,
                p->xsite,
                LCHROM+2,
                cstr,
                p->fitness
               );
    }
}

int main()
{
    SGA_params_t params;
    SGA_stats_t  stats;
    SGA_indiv_t  oldpop[POPSIZE], newpop[POPSIZE];
    FILE *fp1, *fp2;
    int ncross, nmutation;

    fp1 = fopen("gen1.txt", "w");
    fp2 = fopen("gen2.txt", "w");

    params.popsize = POPSIZE;
    params.lchrom  = 30;
    params.pcross  = 0.6;
    params.pmutate = 0.000333;

    SGA_init(123);
    SGA_init_pop(oldpop, &stats, &params);
    SGA_advance(oldpop, newpop, &stats, &params, &ncross, &nmutation);

    report(fp1, oldpop);
    report(fp2, newpop);
}
