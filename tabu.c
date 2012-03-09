#include <stdio.h>
#include <stdlib.h>
#include "tabu.h"

/*
 * computes 2^e mod HASHLEN, where e is the integer represented by the binary
 * string a of length len
 */
int hash(int a[], int len)
{
    int base = 2, result = 1;

    while (--len)
    {
        if (a[len] == 1) result = (result * base) % HASHLEN;
        base = (base * base) % HASHLEN;
    }

    if (a[len] == 1) result = (result * base) % HASHLEN;

    return result;
}

void tabu_init(tabu_t *tl, int size)
{
    tl->size  = size;
    tl->index = 0;
}

void tabu_add(tabu_t *tl, tabu_elem_t *el)
{
    tl->hist[tl->index] = *el;
    tl->index++;
    /* ... */
}

int main()
{
    char filename[256];

    return EXIT_SUCCESS;
}
