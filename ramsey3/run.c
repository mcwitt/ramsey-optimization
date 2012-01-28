/* 
 * File: run.c
 *
 * Version 3
 *
 * Author: Matt Wittmann <mwittman@ucsc.edu>
 *
 * See ramsey.c for description
 *
 * Main simulation loop
 */

#include <limits.h>
#include "ramsey.h"

void run()
{
    rep_t *p;
    char filename[256];
    int it, done;

    min = INT_MAX;
    nsweeps = 0;
    done = 0;
#ifndef NOTIME
    start = clock();
#endif

    while (! done && nsweeps < MAX_SWEEPS)
    {
        sweep();
        nsweeps++;

        temper();

        if (nsweeps % WRITE_INTERVAL == 0)
        {
            sprintf(filename, "%d-%d-%d_%d.bin",
                    S, S, NV, rseed);
            save_state(filename);
            printf("state saved to %s\n", filename);
        }

        for (it = 0; it < nt; it++)
        {
            p = &reps[ri[it]];
            if (p->en < min)
            {
                min = p->en;
                print_status();
                sprintf(filename, "%d-%d-%d_%d.graph", S, S, NV, rseed);
                save_graph(p->sp, filename);

                if (p->en == 0) { done = 1; break; }
            }
        }
    }

    print_status();
}
