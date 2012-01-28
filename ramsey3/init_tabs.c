/* 
 * File: init_tabs.c
 *
 * Version 3
 *
 * Author: Matt Wittmann <mwittman@ucsc.edu>
 *
 * See ramsey.c for description
 *
 * Initializes look-up tables `sub' and `edg'. For a given edge index ei,
 * sub[ei] is an array of the N_v-2 choose S-2 subgraphs that include the edge
 * ei. For a given subgraph si, edg[si] is an array of the S(S-1)/2 edges
 * included in the subgraph si.
 */

#include "ramsey.h"

void init_tabs()
{
    int ps[NED];
    int pe[NSG];
    int c[S+2];     /* array of vertices of the current subgraph */
    int ei, si;     /* edge index, subgraph index */
    int j, k;

    for (j = 0; j < NED; j++)
    {
        sub[j] = (int*) malloc(NSGFE * sizeof(int));
        ps[j] = 0;
    }

    for (j = 0; j < NSG; j++)
    {
        edg[j] = (int*) malloc(neds * sizeof(int));
        pe[j] = 0;
    }

    /* 
     * iterate over all subgraphs with S vertices
     */

    /*
     * algorithm to generate combinations adapted from Algorithm L in Knuth's
     * Art of Computer Programming Vol. 4, Fasc. 3 (all-caps labels
     * correspond to labels in the book)
     */

    /* INITIALIZE */
    si = 0;
    c[S] = NV;
    c[S+1] = 0;
    for (j = 0; j < S; j++) c[j] = j;

    while (1)
    {
        /*
         * VISIT combination c_1 c_2 ... c_S
         * (algorithm guarantees that c_1 < c_2 < ... < c_S)
         */

        /* iterate over edges in this subgraph */
        for (k = 0; k < S; k++)
        {
            for (j = 0; j < k; j++)
            {
                ei = c[k]*(c[k]-1)/2 + c[j];

                /*
                 * add subgraph si to list for edge ei,
                 * and edge ei to list for subgraph si
                 */

                sub[ei][ps[ei]++] = si;
                edg[si][pe[si]++] = ei;
            }
        }

        si++;   /* finished with this subgraph, increment label */

        /* FIND j */
        j = 0;
        while (c[j] + 1 == c[j+1]) { c[j] = j; j++; }

        /* DONE? */
        if (j == S) break;

        c[j]++;
    }
}
