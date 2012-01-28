/* 
 * File: update_fields.c
 *
 * Version 3
 *
 * Author: Matt Wittmann <mwittman@ucsc.edu>
 *
 * See ramsey.c for description
 *
 * Update local fields of all affected spins after a flip is made
 */

#include "ramsey.h"

void update_fields(int ei, int sp[], int nb[], int h2[])
{
    int si, j, ej, nbf;

    if (sp[ei] == 1)
    {
        for (si = 0; si < NSGFE; si++)
        {
            nbf = nb[sub[ei][si]] += 1;
            if (nbf == neds || nbf == 1)
            {
                /* completed a blue clique or destroyed a red clique */
                h2[ei] += 1;
                for (j = 0; j < neds; j++)
                    h2[edg[sub[ei][si]][j]] -= 1;
            }
            else if (nbf == nedsm1)
            {
                /* completed blue clique except for one red edge--
                 * update field at red edge */
                for (j = 0; j < neds; j++)
                {
                    ej = edg[sub[ei][si]][j];
                    if (sp[ej] == -1) { h2[ej] -= 1; break; }
                }
            }
            else if (nbf == 2)
            {
                /* destroyed an almost-complete red clique--
                 * update field at the preexisting blue edge */
                for (j = 0; j < neds; j++)
                {
                    ej = edg[sub[ei][si]][j];
                    if (sp[ej] == 1 && ej != ei) { h2[ej] -= 1; break; }
                }
            }
        }
    }
    else    /* analagous procedure if the new edge is red... */
    {
        for (si = 0; si < NSGFE; si++)
        {
            nbf = nb[sub[ei][si]] -= 1;
            if (nbf == 0 || nbf == nedsm1)
            {
                h2[ei] -= 1;
                for (j = 0; j < neds; j++)
                    h2[edg[sub[ei][si]][j]] += 1;
            }
            else if (nbf == 1)
            {
                for (j = 0; j < neds; j++)
                {
                    ej = edg[sub[ei][si]][j];
                    if (sp[ej] == 1) { h2[ej] += 1; break; }
                }
            }
            else if (nbf == nedsm2)
            {
                for (j = 0; j < neds; j++)
                {
                    ej = edg[sub[ei][si]][j];
                    if (sp[ej] == -1 && ej != ei) { h2[ej] += 1; break; }
                }
            }
        }
    }
}
