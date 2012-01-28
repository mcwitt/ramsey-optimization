/* 
 * File: update_fields.c
 *
 * Version 4
 *
 * Author: Matt Wittmann <mwittman@ucsc.edu>
 *
 * See ramsey.c for description
 *
 * Update local fields of all affected spins after a flip is made
 */

#include "ramsey.h"

void update_fields(int ei, int sp[], int nbr[], int nbs[], int h2[])
{
    int si, j, ej, nbf;

    if (sp[ei] == 1)
    {
        for (si = 0; si < NSGFES; si++)
        {
            nbf = nbs[subs[ei][si]] += 1;

            if (nbf == neds)    /* created blue clique */
            {
                h2[ei] += 1;
                for (j = 0; j < neds; j++) h2[edgs[subs[ei][si]][j]] -= 1;
            }
            else if (nbf == nedsm1) /* created incomplete blue clique */
            {
                for (j = 0; j < neds; j++)
                {
                    ej = edgs[subs[ei][si]][j];
                    if (sp[ej] == -1) { h2[ej] -= 1; break; }
                }
            }
        }

        for (si = 0; si < NSGFER; si++)
        {
            nbf = nbr[subr[ei][si]] += 1;

            if (nbf == 1)   /* destroyed red clique */
            {
                h2[ei] += 1;
                for (j = 0; j < nedr; j++) h2[edgr[subr[ei][si]][j]] -= 1;
            }
            else if (nbf == 2)   /* destroyed incomplete red clique */
            {
                for (j = 0; j < nedr; j++)
                {
                    ej = edgr[subr[ei][si]][j];
                    if (sp[ej] == 1 && ej != ei) { h2[ej] -= 1; break; }
                }
            }
        }
    }
    else
    {
        for (si = 0; si < NSGFER; si++)
        {
            nbf = nbr[subr[ei][si]] -= 1;

            if (nbf == 0)   /* created a red clique */
            {
                h2[ei] -= 1;
                for (j = 0; j < nedr; j++) h2[edgr[subr[ei][si]][j]] += 1;
            }
            else if (nbf == 1)  /* created an incomplete red clique */
            {
                for (j = 0; j < nedr; j++)
                {
                    ej = edgr[subr[ei][si]][j];
                    if (sp[ej] == 1) { h2[ej] += 1; break; }
                }
            }
        }

        for (si = 0; si < NSGFES; si++)
        {
            nbf = nbs[subs[ei][si]] -= 1;

            if (nbf == nedsm1)  /* destroyed a blue clique */
            {
                h2[ei] -= 1;
                for (j = 0; j < neds; j++) h2[edgs[subs[ei][si]][j]] += 1;
            }
            else if (nbf == nedsm2) /* destroyed an incomplete blue clique */
            {
                for (j = 0; j < neds; j++)
                {
                    ej = edgs[subs[ei][si]][j];
                    if (sp[ej] == -1 && ej != ei) { h2[ej] += 1; break; }
                }
            }
        }
    }
}
