#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include "ramsey.h"

#define N_ITER 100

R_replica_t r;
char msg[500];

int check_energy()
{
    int i, R_en;

    for (i = 0; i < N_ITER; i++)
    {
        R_randomize(&r, 0.5, 0);
        R_en = R_energy(r.sp);
        if (r.en != R_en)
        {
            sprintf(msg, "check_energy: Failed on iteration %d. "\
                    "(r.en = %d, R_energy = %d)\n", i, r.en, R_en);
            return 1;
        }
    }

    return 0;
}

int main()
{
    int exit_code;
    uint32_t seed;

    seed = time(NULL);

    R_init(seed);
    R_init_replica(&r);

    exit_code = check_energy();
    if (exit_code != 0) printf("%s\n", msg);
    else printf("checks: all tests successful\n");

    R_finalize();
    return exit_code;
}
