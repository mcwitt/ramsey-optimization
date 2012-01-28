/* 
 * File: temper.c
 *
 * Version 4
 *
 * Author: Matt Wittmann <mwittman@ucsc.edu>
 *
 * See ramsey.c for description
 *
 * Attempt a swap between each pair of temperatures
 */

#include <math.h>
#include "ramsey.h"

void temper()
{
    double logar;
    int it, copy;

    for (it = 1; it < nt; it++)
    {
        logar = (reps[ri[it-1]].en - reps[ri[it]].en)
            * (mbeta[it] - mbeta[it-1]);

        if (URAND() < exp(logar))
        {
            /* do PT swap */
            copy = ri[it-1];
            ri[it-1] = ri[it];
            ri[it] = copy;
            nswaps[it]++;
        }
    }
}
