#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    int nv, ngraph, seed;
    int ned, igraph, j;

    if (argc != 4)
    {
        fprintf(stderr, "Usage: %s nv ngraph seed\n", argv[0]);
        fprintf(stderr, "Writes ngraph random graphs to stdout\n");
        exit(EXIT_FAILURE);
    }

    nv = atoi(argv[1]);
    ngraph = atoi(argv[2]);
    seed = atoi(argv[3]);

    srand(seed);
    ned = nv*(nv-1)/2;

    for (igraph = 0; igraph < ngraph; igraph++)
    {
        printf("%d\n", nv);
        for (j = 0; j < ned; j++) printf("%d\n", rand() % 2);
    }

    return EXIT_SUCCESS;
}
