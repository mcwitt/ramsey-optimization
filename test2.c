#include <stdlib.h>
#include <stdio.h>
#include "qselect.h"

#define MAXLEN QSELECT_MAXLEN

int cmp(const void *a, const void *b) {
    return *(int*)a - *(int*)b;
}

int cmp_index(void *thunk, const void *a, const void *b) {
    int *arr = (int*) thunk;
    int ia = *(int*)a, ib = *(int*)b;
    return arr[ia] - arr[ib];
}

double test_select(int k, int n, int arr[])
{
    qsort(arr, n, sizeof(int), cmp);
    return arr[k];
}

int test_select_index(int k, int n, int arr[])
{
    int i, index[MAXLEN];

    i = n; while (--i) index[i] = i; index[0] = 0;
    qsort_r(index, n, sizeof(int), arr, cmp_index);
    return index[k];
}

int main(int argc, char *argv[])
{
    int i1, i2, r1, r2, a1[MAXLEN], a2[MAXLEN];
    int i, n, k;
    FILE *fp;

    fp = fopen("data.txt", "r");
    
    k = atoi(argv[1]);

    for (i = 0; i < MAXLEN; i++)
        if (fscanf(fp, "%d", &a1[i]) == EOF) break;

    fclose(fp);

    n = i;
    while (--i) a2[i] = a1[i];
    a2[0] = a1[0];

/*    r1 = qselect(k, n, a1);
    r2 = test_select(k, n, a2);
    printf("%d\t%d\n", r1, r2);*/

    i1 = qselect_index(k, n, a1);
    i2 = test_select_index(k, n, a2);
    printf("%d\t%d\n", i1, i2);

    return EXIT_SUCCESS;
}

