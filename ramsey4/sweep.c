/* 
 * File: sweep.c
 *
 * Version 4
 *
 * Author: Matt Wittmann <mwittman@ucsc.edu>
 *
 * See ramsey.c for description
 *
 * Update each spin once using Metropolis probability. When a flip is made,
 * call `update_fields(..)'
 */

#include <math.h>
#include "ramsey.h"

void sweep()
{
    rep_t *p;
    int it, j, delta;

    for (it = 0; it < nt; it++)
    {
        p = &reps[ri[it]];
        
        for (j = 0; j < NED; j++)
        {
            /* compute energy difference of flip */
            delta = p->h2[j]*p->sp[j];

            /* flip with Metropolis probability */
            if (delta <= 0 || URAND() < exp(mbeta[it]*delta))
            {
                p->sp[j] *= -1;
                p->en += delta;
                update_fields(j, p->sp, p->nbr, p->nbs, p->h2);
            }
        }
#ifdef DEBUG
        assert(debug_energy(p->sp) == p->energy);
#endif
    }   /* end of loop over temperatures */
}
