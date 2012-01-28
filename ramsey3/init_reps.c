/* 
 * File: init_reps.c
 *
 * Version 3
 *
 * Author: Matt Wittmann <mwittman@ucsc.edu>
 *
 * See ramsey.c for description
 *
 * Initializes each replica in a random state with equal numbers of red and
 * blue edges on average
 */

#include "ramsey.h"

void init_reps()
{
    rep_t *p;
    int iT, j;

    for (iT = 0; iT < nt; iT++)
    {
        ri[iT] = iT;
        p = &reps[iT];
        p->en = NSG;
        nswaps[iT] = 0;

        /* first initialize in a self-consistent state with all edges blue */
        for (j = 0; j < NSG; j++) p->nb[j] = neds;

        for (j = 0; j < NED; j++)
        {
            p->sp[j] = 1;
            p->h2[j] = -NSGFE;
        }

        /* randomize spins, updating fields with each flip */
        for (j = 0; j < NED; j++)
        {
            if (URAND() < 0.5)
            {
                p->sp[j] = -1;
                p->en += p->h2[j];
                update_fields(j, p->sp, p->nb, p->h2);
            }
        }
    }
}
