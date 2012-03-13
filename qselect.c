/* Code adapted from Numerical Recipes in C (2nd & 3rd Eds.) */

#include "qselect.h"

#define SWAP(a,b) { temp=(a); (a)=(b); (b)=temp; }

double qselect(unsigned long k, unsigned long n, double arr[])
/*
 * Given k in [0..n-1] returns an array value from arr[0..n-1] such that k
 * array values are less than or equal to the one returned. The input array
 * will be rearranged to have this value in location arr[k], with all smaller
 * elements movd to arr[0..k-1] (in arbitrary order) and all larger elements in
 * arr[k+1..n-1] (also in arbitrary order).
 */
{
    unsigned long i, ir, j, l, mid;
    double a, temp;

    l = 0;
    ir = n-1;

    for (;;) {
        if (ir <= l+1) {    /* Active partition contains 1 or 2 elements */
            if (ir == l+1 && arr[ir] < arr[l])  /* Case of 2 elements */
                SWAP(arr[l], arr[ir]);
            return arr[k];
        } else {
            /*
             * Choose median of left, center, and right elements as
             * partitioning element a. Also rearrange so that
             * arr[l] <= arr[l+1], arr[ir] >= arr[l+1].
             */
            mid = (l+ir) >> 1;
            SWAP(arr[mid], arr[l+1]);
            if (arr[l] > arr[ir])
                SWAP(arr[l], arr[ir]);
            if (arr[l+1] > arr[ir])
                SWAP(arr[l+1], arr[ir]);
            if (arr[l] > arr[l+1])
                SWAP(arr[l], arr[l+1]);

            i = l+1;                /* Initialize pointers for partitioning */
            j = ir;
            a = arr[l+1];           /* Partitioning element */
            for (;;) {              /* Beginning of innermost loop */
                do i++; while (arr[i] < a); /* Scan up to find element > a */
                do j--; while (arr[j] > a); /* Scan down to find element < a */
                if (j < i) break;   /* Pointers crossed. Done partitioning */
                SWAP(arr[i], arr[j]);
            }                       /* End of innermost loop */
            arr[l+1] = arr[j];      /* Insert partitioning element */
            arr[j] = a;
            if (j >= k) ir = j-1;   /* Keep active the partition that contains */
            if (j <= k) l = i;      /*   the kth element */
        }
    }
}

unsigned long qselect_index(unsigned long k, unsigned long n, double arr[])
/* Returns the index of the kth smallest value in the array arr[1..n] */
{
    unsigned long i, ir, j, l, mid, temp, tab[QSELECT_MAXLEN];
    double a;

    l = 0;
    ir = n-1;

    while (--n) tab[n] = n;   /* Initialize index table */
    tab[0] = 0;

    for (;;) {
        if (ir <= l+1) {    /* Active partition contains 1 or 2 elements */
            if (ir == l+1 && arr[tab[ir]] < arr[tab[l]]) /* Case of 2 */
                SWAP(tab[l], tab[ir]);
            return tab[k];
        } else {
            /*
             * Choose median of left, center, and right elements as
             * partitioning element a. Also rearrange so that arr[l] <=
             * arr[l+1], arr[ir] >= arr[l+1].
             */
            mid = (l+ir) >> 1;
            SWAP(tab[mid], tab[l+1]);
            if (arr[tab[l]] > arr[tab[ir]])
                SWAP(tab[l], tab[ir]);
            if (arr[tab[l+1]] > arr[tab[ir]])
                SWAP(tab[l+1], tab[ir]);
            if (arr[tab[l]] > arr[tab[l+1]])
                SWAP(tab[l], tab[l+1]);

            i = l+1;                /* Initialize pointers for partitioning */
            j = ir;
            a = arr[tab[l+1]];      /* Partitioning element */
            for (;;) {              /* Beginning of innermost loop */
                do i++; while (arr[tab[i]] < a); /* Scan up to find # > a */
                do j--; while (arr[tab[j]] > a); /* Scan down to find # < a */
                if (j < i) break;   /* Pointers crossed. Done partitioning */
                SWAP(tab[i], tab[j]);
            }                       /* End of innermost loop */
            SWAP(tab[l+1], tab[j]); /* Insert partitioning element */
            if (j >= k) ir = j-1;   /* Keep active the partition that contains */
            if (j <= k) l = i;      /*   the kth element */
        }
    }
}
