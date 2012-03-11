/* Code adapted from Numerical Recipes in C (2nd Ed.) */

#define SWAP(a,b) { temp=(a); (a)=(b); (b)=temp; }

double select(unsigned long k, unsigned long n, double arr[])
/*
 * Returns the kth smallest value in the array arr[1..n]. The input array will
 * be rearranged to have this value in location arr[k], with all smaller
 * elements moved to arr[1..k-1] (in arbitrary order) and all larger elements
 * in arr[k+1..n] (also in arbitrary order).
 */
{
    unsigned long i, ir, j, l, mid;
    double a, temp;

    l = 1;
    ir = n;

    for (;;) {
        if (ir <= l+1) {    /* Active partition contains 1 or 2 elements */
            if (ir == l+1 && arr[ir] < arr[l])  /* Case of 2 elements */
                SWAP(arr[l], arr[ir]);
            return arr[k];
        } else {
            /*
             * Choose median of left, center, and right elements as
             * partitioning element a. Also rearrange so that arr[l] <=
             * arr[l+1], arr[ir] >= arr[l+1].
             */
            mid = (l+ir) >> 1;
            SWAP(arr[mid], arr[l+1]);
            if (arr[l  ] > arr[ir ]) SWAP(arr[l  ], arr[ir ]);
            if (arr[l+1] > arr[ir ]) SWAP(arr[l+1], arr[ir ]);
            if (arr[l  ] > arr[l+1]) SWAP(arr[l  ], arr[l+1]);
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




