/* 
 * File: io.c
 *
 * Version 4
 *
 * Author: Matt Wittmann <mwittman@ucsc.edu>
 *
 * See ramsey.c for description
 *
 * Provides file I/O operations.
 */

#include "ramsey.h"

void save_graph(int sp[NED], char filename[])
{
    FILE *fp;
    int i;

    fp = fopen(filename, "w");
    fprintf(fp, "%d\n", NV);
    fprintf(fp, "%d\n", R);
    fprintf(fp, "%d\n", S);
    
    for (i = 0; i < NED; i++)
        fprintf(fp, "%d\n", (sp[i] == 1) ? 1 : 0);

    fclose(fp);
}

void save_state(char filename[])
{
    FILE *fp;

    fp = fopen(filename, "w");
    fwrite(ri, sizeof(int), nt, fp);
    fwrite(reps, sizeof(rep_t), nt, fp);
    fclose(fp);
}

void load_state(char filename[])
{
    FILE *fp;

    fp = fopen(filename, "r");
    assert(fread(ri, sizeof(int), nt, fp) == nt);
    assert(fread(reps, sizeof(rep_t), nt, fp) == nt);
    fclose(fp);
}

void print_status()
{
    int it;
#ifndef NOTIME
    int trun;
#endif

    printf("\n");
    printf("min. energy     : %d\n", min);
    printf("# of sweeps     : %d\n", nsweeps);
#ifndef NOTIME
    trun = (clock() - start)/CLOCKS_PER_SEC;
    printf("time running    : %d seconds\n", trun);
    printf("sweep rate      : %.2f / s\n", (float) nsweeps/trun);
#endif
    for (it = 0; it < nt; it++)
        printf("%5d ", reps[ri[it]].en);
    printf("\n");
    for (it = 0; it < nt; it++)
        printf("%5.2f ", T[it]);
    printf("\n   ");
    for (it = 1; it < nt; it++)
        printf("%5.2f ", (float) nswaps[it]/nsweeps);
    printf("\n");
    fflush(stdout);
}
