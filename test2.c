#include <stdlib.h>
#include <stdio.h>
#include "quickselect.h"

#define MAXLEN 2000

int cmp(const void *a, const void *b) {
    return (*(double*)a < *(double*)b) ? -1: 1;
}

double test_select(int k, int n, double arr[])
{
    qsort(&arr[1], n, sizeof(double), cmp);
    return arr[k];
}

int test_select_index(int k, int n, double arr[])
{
}

int main(int argc, char *argv[])
{
    double r1, r2, a1[MAXLEN], a2[MAXLEN];
    int i, n, k;
    FILE *fp;

    fp = fopen("data.txt", "r");
    
    k = atoi(argv[1]);

    for (i = 1; i < MAXLEN; i++)
        if (fscanf(fp, "%lf", &a1[i]) == EOF) break;

    fclose(fp);

    n = i-1;
    while (--i) a2[i] = a1[i];

    r1 = quickselect(k, n, a1);
    r2 = test_select(k, n, a2);

    printf("%f\t%f\n", r1, r2);

    return EXIT_SUCCESS;
}

